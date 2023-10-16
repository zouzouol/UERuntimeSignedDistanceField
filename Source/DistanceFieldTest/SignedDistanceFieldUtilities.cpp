// Fill out your copyright notice in the Description page of Project Settings.


#include "SignedDistanceFieldUtilities.h"

#include "Kismet/GameplayStatics.h"
#include "Engine/StaticMeshActor.h"
#include "DistanceFieldAtlas.h"
#include "MeshCardRepresentation.h"
#include "ObjectCacheContext.h"

static FVector3f UniformSampleHemisphere(FVector2D Uniforms)
{
	Uniforms = Uniforms * 2.0f - 1.0f;

	if (Uniforms == FVector2D::ZeroVector)
	{
		return FVector3f::ZeroVector;
	}

	float R;
	float Theta;

	if (FMath::Abs(Uniforms.X) > FMath::Abs(Uniforms.Y))
	{
		R = Uniforms.X;
		Theta = (float)PI / 4 * (Uniforms.Y / Uniforms.X);
	}
	else
	{
		R = Uniforms.Y;
		Theta = (float)PI / 2 - (float)PI / 4 * (Uniforms.X / Uniforms.Y);
	}

	// concentric disk sample
	const float U = R * FMath::Cos(Theta);
	const float V = R * FMath::Sin(Theta);
	const float R2 = R * R;

	// map to hemisphere [P. Shirley, Kenneth Chiu; 1997; A Low Distortion Map Between Disk and Square]
	return FVector3f(U * FMath::Sqrt(2 - R2), V * FMath::Sqrt(2 - R2), 1.0f - R2);
}

#if USE_EMBREE
void EmbreeFilterFunc(const struct RTCFilterFunctionNArguments* args)
{
	FEmbreeGeometry* EmbreeGeometry = (FEmbreeGeometry*)args->geometryUserPtr;
	FEmbreeTriangleDesc Desc = EmbreeGeometry->TriangleDescs[RTCHitN_primID(args->hit, 1, 0)];

	FEmbreeIntersectionContext& IntersectionContext = *static_cast<FEmbreeIntersectionContext*>(args->context);
	IntersectionContext.ElementIndex = Desc.ElementIndex;

	const RTCHit& EmbreeHit = *(RTCHit*)args->hit;
	if (IntersectionContext.SkipPrimId != RTC_INVALID_GEOMETRY_ID && IntersectionContext.SkipPrimId == EmbreeHit.primID)
	{
		// Ignore hit in order to continue tracing
		args->valid[0] = 0;
	}
}

void EmbreeErrorFunc(void* userPtr, RTCError code, const char* str)
{
	FString ErrorString;
	TArray<TCHAR, FString::AllocatorType>& ErrorStringArray = ErrorString.GetCharArray();
	ErrorStringArray.Empty();

	int32 StrLen = FCStringAnsi::Strlen(str);
	int32 Length = FUTF8ToTCHAR_Convert::ConvertedLength(str, StrLen);
	ErrorStringArray.AddUninitialized(Length + 1); // +1 for the null terminator
	FUTF8ToTCHAR_Convert::Convert(ErrorStringArray.GetData(), ErrorStringArray.Num(), reinterpret_cast<const ANSICHAR*>(str), StrLen);
	ErrorStringArray[Length] = TEXT('\0');

	UE_LOG(LogTemp, Error, TEXT("Embree error: %s Code=%u"), *ErrorString, (uint32)code);
}
#endif

USignedDistanceFieldUtilities::USignedDistanceFieldUtilities()
{
}

void USignedDistanceFieldUtilities::GenerateStratifiedUniformHemisphereSamples(int32 NumSamples, FRandomStream& RandomStream, TArray<FVector3f>& Samples)
{
	const int32 NumSamplesDim = FMath::TruncToInt(FMath::Sqrt((float)NumSamples));

	Samples.Empty(NumSamplesDim * NumSamplesDim);

	for (int32 IndexX = 0; IndexX < NumSamplesDim; IndexX++)
	{
		for (int32 IndexY = 0; IndexY < NumSamplesDim; IndexY++)
		{
			const float U1 = RandomStream.GetFraction();
			const float U2 = RandomStream.GetFraction();

			const float Fraction1 = (IndexX + U1) / (float)NumSamplesDim;
			const float Fraction2 = (IndexY + U2) / (float)NumSamplesDim;

			Samples.Add(UniformSampleHemisphere(FVector2D(Fraction1, Fraction2)));
		}
	}
}

FMatrix44f USignedDistanceFieldUtilities::GetTangentBasisFrisvad(FVector3f TangentZ)
{
	FVector3f TangentX;
	FVector3f TangentY;

	if (TangentZ.Z < -0.9999999f)
	{
		TangentX = FVector3f(0, -1, 0);
		TangentY = FVector3f(-1, 0, 0);
	}
	else
	{
		float A = 1.0f / (1.0f + TangentZ.Z);
		float B = -TangentZ.X * TangentZ.Y * A;
		TangentX = FVector3f(1.0f - TangentZ.X * TangentZ.X * A, B, -TangentZ.X);
		TangentY = FVector3f(B, 1.0f - TangentZ.Y * TangentZ.Y * A, -TangentZ.Y);
	}

	FMatrix44f LocalBasis;
	LocalBasis.SetIdentity();
	LocalBasis.SetAxis(0, TangentX);
	LocalBasis.SetAxis(1, TangentY);
	LocalBasis.SetAxis(2, TangentZ);
	return LocalBasis;
}

void USignedDistanceFieldUtilities::SetupEmbreeScene(FString MeshName, const FSourceMeshDataForDerivedDataTask& SourceMeshData, const FStaticMeshLODResources& LODModel, const TArray<FSignedDistanceFieldBuildMaterialData>& MaterialBlendModes, bool bGenerateAsIfTwoSided, FEmbreeScene& EmbreeScene)
{
	const uint32 NumIndices = SourceMeshData.IsValid() ? SourceMeshData.GetNumIndices() : LODModel.IndexBuffer.GetNumIndices();
	const int32 NumTriangles = NumIndices / 3;
	const uint32 NumVertices = SourceMeshData.IsValid() ? SourceMeshData.GetNumVertices() : LODModel.VertexBuffers.PositionVertexBuffer.GetNumVertices();
	EmbreeScene.NumIndices = NumTriangles;

	TArray<FkDOPBuildCollisionTriangle<uint32> > BuildTriangles;

#if USE_EMBREE
	EmbreeScene.bUseEmbree = true;

	if (EmbreeScene.bUseEmbree)
	{
		EmbreeScene.EmbreeDevice = rtcNewDevice(nullptr);
		rtcSetDeviceErrorFunction(EmbreeScene.EmbreeDevice, EmbreeErrorFunc, nullptr);

		RTCError ReturnErrorNewDevice = rtcGetDeviceError(EmbreeScene.EmbreeDevice);
		if (ReturnErrorNewDevice != RTC_ERROR_NONE)
		{
			UE_LOG(LogTemp, Warning, TEXT("GenerateSignedDistanceFieldVolumeData failed for %s. Embree rtcNewDevice failed. Code: %d"), *MeshName, (int32)ReturnErrorNewDevice);
			return;
		}

		EmbreeScene.EmbreeScene = rtcNewScene(EmbreeScene.EmbreeDevice);
		rtcSetSceneFlags(EmbreeScene.EmbreeScene, RTC_SCENE_FLAG_NONE);

		RTCError ReturnErrorNewScene = rtcGetDeviceError(EmbreeScene.EmbreeDevice);
		if (ReturnErrorNewScene != RTC_ERROR_NONE)
		{
			UE_LOG(LogTemp, Warning, TEXT("GenerateSignedDistanceFieldVolumeData failed for %s. Embree rtcNewScene failed. Code: %d"), *MeshName, (int32)ReturnErrorNewScene);
			rtcReleaseDevice(EmbreeScene.EmbreeDevice);
			return;
		}
	}
#endif

	TArray<int32> FilteredTriangles;
	FilteredTriangles.Empty(NumTriangles);

	if (SourceMeshData.IsValid())
	{
		for (int32 TriangleIndex = 0; TriangleIndex < NumTriangles; ++TriangleIndex)
		{
			const uint32 I0 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 0];
			const uint32 I1 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 1];
			const uint32 I2 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 2];

			const FVector3f V0 = SourceMeshData.VertexPositions[I0];
			const FVector3f V1 = SourceMeshData.VertexPositions[I1];
			const FVector3f V2 = SourceMeshData.VertexPositions[I2];

			const FVector3f TriangleNormal = ((V1 - V2) ^ (V0 - V2));
			const bool bDegenerateTriangle = TriangleNormal.SizeSquared() < SMALL_NUMBER;
			if (!bDegenerateTriangle)
			{
				FilteredTriangles.Add(TriangleIndex);
			}
		}
	}
	else
	{
		for (int32 TriangleIndex = 0; TriangleIndex < NumTriangles; ++TriangleIndex)
		{
			const FIndexArrayView Indices = LODModel.IndexBuffer.GetArrayView();
			const uint32 I0 = Indices[TriangleIndex * 3 + 0];
			const uint32 I1 = Indices[TriangleIndex * 3 + 1];
			const uint32 I2 = Indices[TriangleIndex * 3 + 2];

			const FVector3f V0 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I0);
			const FVector3f V1 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I1);
			const FVector3f V2 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I2);

			const FVector3f TriangleNormal = ((V1 - V2) ^ (V0 - V2));
			const bool bDegenerateTriangle = TriangleNormal.SizeSquared() < SMALL_NUMBER;
			if (!bDegenerateTriangle)
			{
				bool bTriangleIsOpaqueOrMasked = false;

				for (int32 SectionIndex = 0; SectionIndex < LODModel.Sections.Num(); SectionIndex++)
				{
					const FStaticMeshSection& Section = LODModel.Sections[SectionIndex];

					if ((uint32)(TriangleIndex * 3) >= Section.FirstIndex && (uint32)(TriangleIndex * 3) < Section.FirstIndex + Section.NumTriangles * 3)
					{
						if (MaterialBlendModes.IsValidIndex(Section.MaterialIndex))
						{
							bTriangleIsOpaqueOrMasked = !IsTranslucentBlendMode(MaterialBlendModes[Section.MaterialIndex].BlendMode) && MaterialBlendModes[Section.MaterialIndex].bAffectDistanceFieldLighting;
						}

						break;
					}
				}

				if (bTriangleIsOpaqueOrMasked)
				{
					FilteredTriangles.Add(TriangleIndex);
				}
			}
		}
	}

	const int32 NumBufferVerts = 1; // Reserve extra space at the end of the array, as embree has an internal bug where they read and discard 4 bytes off the end of the array
	EmbreeScene.Geometry.VertexArray.Empty(NumVertices + NumBufferVerts);
	EmbreeScene.Geometry.VertexArray.AddUninitialized(NumVertices + NumBufferVerts);

	const int32 NumFilteredIndices = FilteredTriangles.Num() * 3;

	EmbreeScene.Geometry.IndexArray.Empty(NumFilteredIndices);
	EmbreeScene.Geometry.IndexArray.AddUninitialized(NumFilteredIndices);

	FVector3f* EmbreeVertices = EmbreeScene.Geometry.VertexArray.GetData();
	uint32* EmbreeIndices = EmbreeScene.Geometry.IndexArray.GetData();
	EmbreeScene.Geometry.TriangleDescs.Empty(FilteredTriangles.Num());

	for (int32 FilteredTriangleIndex = 0; FilteredTriangleIndex < FilteredTriangles.Num(); FilteredTriangleIndex++)
	{
		uint32 I0, I1, I2;
		FVector3f V0, V1, V2;

		const int32 TriangleIndex = FilteredTriangles[FilteredTriangleIndex];
		if (SourceMeshData.IsValid())
		{
			I0 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 0];
			I1 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 1];
			I2 = SourceMeshData.TriangleIndices[TriangleIndex * 3 + 2];

			V0 = SourceMeshData.VertexPositions[I0];
			V1 = SourceMeshData.VertexPositions[I1];
			V2 = SourceMeshData.VertexPositions[I2];
		}
		else
		{
			const FIndexArrayView Indices = LODModel.IndexBuffer.GetArrayView();
			I0 = Indices[TriangleIndex * 3 + 0];
			I1 = Indices[TriangleIndex * 3 + 1];
			I2 = Indices[TriangleIndex * 3 + 2];

			V0 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I0);
			V1 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I1);
			V2 = LODModel.VertexBuffers.PositionVertexBuffer.VertexPosition(I2);
		}

		bool bTriangleIsTwoSided = false;

		for (int32 SectionIndex = 0; SectionIndex < LODModel.Sections.Num(); SectionIndex++)
		{
			const FStaticMeshSection& Section = LODModel.Sections[SectionIndex];

			if ((uint32)(TriangleIndex * 3) >= Section.FirstIndex && (uint32)(TriangleIndex * 3) < Section.FirstIndex + Section.NumTriangles * 3)
			{
				if (MaterialBlendModes.IsValidIndex(Section.MaterialIndex))
				{
					bTriangleIsTwoSided = MaterialBlendModes[Section.MaterialIndex].bTwoSided;
				}

				break;
			}
		}

		if (EmbreeScene.bUseEmbree)
		{
			EmbreeIndices[FilteredTriangleIndex * 3 + 0] = I0;
			EmbreeIndices[FilteredTriangleIndex * 3 + 1] = I1;
			EmbreeIndices[FilteredTriangleIndex * 3 + 2] = I2;

			EmbreeVertices[I0] = V0;
			EmbreeVertices[I1] = V1;
			EmbreeVertices[I2] = V2;

			FEmbreeTriangleDesc Desc;
			// Store bGenerateAsIfTwoSided in material index
			Desc.ElementIndex = bGenerateAsIfTwoSided || bTriangleIsTwoSided ? 1 : 0;
			EmbreeScene.Geometry.TriangleDescs.Add(Desc);
		}
		else
		{
			BuildTriangles.Add(FkDOPBuildCollisionTriangle<uint32>(
				// Store bGenerateAsIfTwoSided in material index
				bGenerateAsIfTwoSided || bTriangleIsTwoSided ? 1 : 0,
				FVector(V0),
				FVector(V1),
				FVector(V2)));
		}
	}

#if USE_EMBREE
	if (EmbreeScene.bUseEmbree)
	{
		RTCGeometry Geometry = rtcNewGeometry(EmbreeScene.EmbreeDevice, RTC_GEOMETRY_TYPE_TRIANGLE);
		EmbreeScene.Geometry.InternalGeometry = Geometry;

		rtcSetSharedGeometryBuffer(Geometry, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, EmbreeVertices, 0, sizeof(FVector3f), NumVertices);
		rtcSetSharedGeometryBuffer(Geometry, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, EmbreeIndices, 0, sizeof(uint32) * 3, FilteredTriangles.Num());

		rtcSetGeometryUserData(Geometry, &EmbreeScene.Geometry);
		rtcSetGeometryIntersectFilterFunction(Geometry, EmbreeFilterFunc);

		rtcCommitGeometry(Geometry);
		rtcAttachGeometry(EmbreeScene.EmbreeScene, Geometry);
		rtcReleaseGeometry(Geometry);

		rtcCommitScene(EmbreeScene.EmbreeScene);

		RTCError ReturnError = rtcGetDeviceError(EmbreeScene.EmbreeDevice);
		if (ReturnError != RTC_ERROR_NONE)
		{
			UE_LOG(LogTemp, Warning, TEXT("GenerateSignedDistanceFieldVolumeData failed for %s. Embree rtcCommitScene failed. Code: %d"), *MeshName, (int32)ReturnError);
			return;
		}
	}
	else
#endif
	{
		EmbreeScene.kDopTree.Build(BuildTriangles);
	}

	// bMostlyTwoSided
	{
		uint32 NumTrianglesTotal = 0;
		uint32 NumTwoSidedTriangles = 0;

		for (int32 SectionIndex = 0; SectionIndex < LODModel.Sections.Num(); SectionIndex++)
		{
			const FStaticMeshSection& Section = LODModel.Sections[SectionIndex];

			if (MaterialBlendModes.IsValidIndex(Section.MaterialIndex))
			{
				NumTrianglesTotal += Section.NumTriangles;

				if (MaterialBlendModes[Section.MaterialIndex].bTwoSided)
				{
					NumTwoSidedTriangles += Section.NumTriangles;
				}
			}
		}

		EmbreeScene.bMostlyTwoSided = NumTwoSidedTriangles * 4 >= NumTrianglesTotal || bGenerateAsIfTwoSided;
	}
}

void USignedDistanceFieldUtilities::DeleteEmbreeScene(FEmbreeScene& EmbreeScene)
{
#if USE_EMBREE
	if (EmbreeScene.bUseEmbree)
	{
		rtcReleaseScene(EmbreeScene.EmbreeScene);
		rtcReleaseDevice(EmbreeScene.EmbreeDevice);
	}
#endif
}





#if USE_EMBREE


class FEmbreePointQueryContext : public RTCPointQueryContext
{
public:
	RTCGeometry MeshGeometry;
	int32 NumTriangles;
};

bool EmbreePointQueryFunction(RTCPointQueryFunctionArguments* args)
{
	const FEmbreePointQueryContext* Context = (const FEmbreePointQueryContext*)args->context;

	check(args->userPtr);
	float& ClosestDistanceSq = *(float*)(args->userPtr);

	const int32 TriangleIndex = args->primID;
	check(TriangleIndex < Context->NumTriangles);

	const FVector3f* VertexBuffer = (const FVector3f*)rtcGetGeometryBufferData(Context->MeshGeometry, RTC_BUFFER_TYPE_VERTEX, 0);
	const uint32* IndexBuffer = (const uint32*)rtcGetGeometryBufferData(Context->MeshGeometry, RTC_BUFFER_TYPE_INDEX, 0);

	const uint32 I0 = IndexBuffer[TriangleIndex * 3 + 0];
	const uint32 I1 = IndexBuffer[TriangleIndex * 3 + 1];
	const uint32 I2 = IndexBuffer[TriangleIndex * 3 + 2];

	const FVector3f V0 = VertexBuffer[I0];
	const FVector3f V1 = VertexBuffer[I1];
	const FVector3f V2 = VertexBuffer[I2];

	const FVector3f QueryPosition(args->query->x, args->query->y, args->query->z);
	const FVector3f ClosestPoint = (FVector3f)FMath::ClosestPointOnTriangleToPoint((FVector)QueryPosition, (FVector)V0, (FVector)V1, (FVector)V2);
	const float QueryDistanceSq = (ClosestPoint - QueryPosition).SizeSquared();

	if (QueryDistanceSq < ClosestDistanceSq)
	{
		ClosestDistanceSq = QueryDistanceSq;

		bool bShrinkQuery = true;

		if (bShrinkQuery)
		{
			args->query->radius = FMath::Sqrt(ClosestDistanceSq);
			// Return true to indicate that the query radius has shrunk
			return true;
		}
	}

	// Return false to indicate that the query radius hasn't changed
	return false;
}

static int32 ComputeLinearVoxelIndex(FIntVector VoxelCoordinate, FIntVector VolumeDimensions)
{
	return (VoxelCoordinate.Z * VolumeDimensions.Y + VoxelCoordinate.Y) * VolumeDimensions.X + VoxelCoordinate.X;
}


class FSparseMeshDistanceFieldAsyncTask
{
public:
	FSparseMeshDistanceFieldAsyncTask(
		const FEmbreeScene& InEmbreeScene,
		const TArray<FVector3f>* InSampleDirections,
		float InLocalSpaceTraceDistance,
		FBox InVolumeBounds,
		float InLocalToVolumeScale,
		FVector2D InDistanceFieldToVolumeScaleBias,
		FIntVector InBrickCoordinate,
		FIntVector InIndirectionSize,
		bool bInUsePointQuery)
		:
		EmbreeScene(InEmbreeScene),
		SampleDirections(InSampleDirections),
		LocalSpaceTraceDistance(InLocalSpaceTraceDistance),
		VolumeBounds(InVolumeBounds),
		LocalToVolumeScale(InLocalToVolumeScale),
		DistanceFieldToVolumeScaleBias(InDistanceFieldToVolumeScaleBias),
		BrickCoordinate(InBrickCoordinate),
		IndirectionSize(InIndirectionSize),
		bUsePointQuery(bInUsePointQuery),
		BrickMaxDistance(MIN_uint8),
		BrickMinDistance(MAX_uint8)
	{}

	void DoWork();

	// Readonly inputs
	const FEmbreeScene& EmbreeScene;
	const TArray<FVector3f>* SampleDirections;
	float LocalSpaceTraceDistance;
	FBox VolumeBounds;
	float LocalToVolumeScale;
	FVector2D DistanceFieldToVolumeScaleBias;
	FIntVector BrickCoordinate;
	FIntVector IndirectionSize;
	bool bUsePointQuery;

	// Output
	uint8 BrickMaxDistance;
	uint8 BrickMinDistance;
	TArray<uint8> DistanceFieldVolume;
};

int32 DebugX = 0;
int32 DebugY = 0;
int32 DebugZ = 0;

void FSparseMeshDistanceFieldAsyncTask::DoWork()
{
	TRACE_CPUPROFILER_EVENT_SCOPE(FSparseMeshDistanceFieldAsyncTask::DoWork);

	const FVector IndirectionVoxelSize = VolumeBounds.GetSize() / FVector(IndirectionSize);
	const FVector DistanceFieldVoxelSize = IndirectionVoxelSize / FVector(DistanceField::UniqueDataBrickSize);
	const FVector BrickMinPosition = VolumeBounds.Min + FVector(BrickCoordinate) * IndirectionVoxelSize;

	DistanceFieldVolume.Empty(DistanceField::BrickSize * DistanceField::BrickSize * DistanceField::BrickSize);
	DistanceFieldVolume.AddZeroed(DistanceField::BrickSize * DistanceField::BrickSize * DistanceField::BrickSize);

	for (int32 ZIndex = 0; ZIndex < DistanceField::BrickSize; ZIndex++)
	{
		for (int32 YIndex = 0; YIndex < DistanceField::BrickSize; YIndex++)
		{
			for (int32 XIndex = 0; XIndex < DistanceField::BrickSize; XIndex++)
			{
				if (XIndex == DebugX && YIndex == DebugY && ZIndex == DebugZ)
				{
					int32 DebugBreak = 0;
				}

				const FVector VoxelPosition = FVector(XIndex, YIndex, ZIndex) * DistanceFieldVoxelSize + BrickMinPosition;
				const int32 Index = (ZIndex * DistanceField::BrickSize * DistanceField::BrickSize + YIndex * DistanceField::BrickSize + XIndex);

				float MinLocalSpaceDistance = LocalSpaceTraceDistance;

				bool bTraceRays = true;

				if (bUsePointQuery)
				{
					RTCPointQuery PointQuery;
					PointQuery.x = VoxelPosition.X;
					PointQuery.y = VoxelPosition.Y;
					PointQuery.z = VoxelPosition.Z;
					PointQuery.time = 0;
					PointQuery.radius = LocalSpaceTraceDistance;

					FEmbreePointQueryContext QueryContext;
					rtcInitPointQueryContext(&QueryContext);
					QueryContext.MeshGeometry = EmbreeScene.Geometry.InternalGeometry;
					QueryContext.NumTriangles = EmbreeScene.Geometry.TriangleDescs.Num();
					float ClosestUnsignedDistanceSq = (LocalSpaceTraceDistance * 2.0f) * (LocalSpaceTraceDistance * 2.0f);
					rtcPointQuery(EmbreeScene.EmbreeScene, &PointQuery, &QueryContext, EmbreePointQueryFunction, &ClosestUnsignedDistanceSq);

					const float ClosestDistance = FMath::Sqrt(ClosestUnsignedDistanceSq);
					bTraceRays = ClosestDistance <= LocalSpaceTraceDistance;
					MinLocalSpaceDistance = FMath::Min(MinLocalSpaceDistance, ClosestDistance);
				}

				if (bTraceRays)
				{
					int32 Hit = 0;
					int32 HitBack = 0;

					for (int32 SampleIndex = 0; SampleIndex < SampleDirections->Num(); SampleIndex++)
					{
						const FVector UnitRayDirection = (FVector)(*SampleDirections)[SampleIndex];
						const float PullbackEpsilon = 1.e-4f;
						// Pull back the starting position slightly to make sure we hit a triangle that VoxelPosition is exactly on.  
						// This happens a lot with boxes, since we trace from voxel corners.
						const FVector StartPosition = VoxelPosition - PullbackEpsilon * LocalSpaceTraceDistance * UnitRayDirection;
						const FVector EndPosition = VoxelPosition + UnitRayDirection * LocalSpaceTraceDistance;

						if (FMath::LineBoxIntersection(VolumeBounds, VoxelPosition, EndPosition, UnitRayDirection))
						{
							FEmbreeRay EmbreeRay;

							FVector RayDirection = EndPosition - VoxelPosition;
							EmbreeRay.ray.org_x = StartPosition.X;
							EmbreeRay.ray.org_y = StartPosition.Y;
							EmbreeRay.ray.org_z = StartPosition.Z;
							EmbreeRay.ray.dir_x = RayDirection.X;
							EmbreeRay.ray.dir_y = RayDirection.Y;
							EmbreeRay.ray.dir_z = RayDirection.Z;
							EmbreeRay.ray.tnear = 0;
							EmbreeRay.ray.tfar = 1.0f;

							FEmbreeIntersectionContext EmbreeContext;
							rtcInitIntersectContext(&EmbreeContext);
							rtcIntersect1(EmbreeScene.EmbreeScene, &EmbreeContext, &EmbreeRay);

							if (EmbreeRay.hit.geomID != RTC_INVALID_GEOMETRY_ID && EmbreeRay.hit.primID != RTC_INVALID_GEOMETRY_ID)
							{
								check(EmbreeContext.ElementIndex != -1);
								Hit++;

								const FVector HitNormal = (FVector)EmbreeRay.GetHitNormal();

								if (FVector::DotProduct(UnitRayDirection, HitNormal) > 0 && !EmbreeContext.IsHitTwoSided())
								{
									HitBack++;
								}

								if (!bUsePointQuery)
								{
									const float CurrentDistance = EmbreeRay.ray.tfar * LocalSpaceTraceDistance;

									if (CurrentDistance < MinLocalSpaceDistance)
									{
										MinLocalSpaceDistance = CurrentDistance;
									}
								}
							}
						}
					}

					// Consider this voxel 'inside' an object if we hit a significant number of backfaces
					if (Hit > 0 && HitBack > .25f * SampleDirections->Num())
					{
						MinLocalSpaceDistance *= -1;
					}
				}

				// Transform to the tracing shader's Volume space
				const float VolumeSpaceDistance = MinLocalSpaceDistance * LocalToVolumeScale;
				// Transform to the Distance Field texture's space
				const float RescaledDistance = (VolumeSpaceDistance - DistanceFieldToVolumeScaleBias.Y) / DistanceFieldToVolumeScaleBias.X;
				check(DistanceField::DistanceFieldFormat == PF_G8);
				const uint8 QuantizedDistance = FMath::Clamp<int32>(FMath::FloorToInt(RescaledDistance * 255.0f + .5f), 0, 255);
				DistanceFieldVolume[Index] = QuantizedDistance;
				BrickMaxDistance = FMath::Max(BrickMaxDistance, QuantizedDistance);
				BrickMinDistance = FMath::Min(BrickMinDistance, QuantizedDistance);
			}
		}
	}
}

void USignedDistanceFieldUtilities::GenerateSignedDistanceFieldVolumeData(
	FString MeshName,
	const FSourceMeshDataForDerivedDataTask& SourceMeshData,
	const FStaticMeshLODResources& LODModel,
	const TArray<FSignedDistanceFieldBuildMaterialData>& MaterialBlendModes,
	const FBoxSphereBounds& Bounds,
	float DistanceFieldResolutionScale,
	bool bGenerateAsIfTwoSided,
	FDistanceFieldVolumeData& OutData)
{
	TRACE_CPUPROFILER_EVENT_SCOPE(GenerateSignedDistanceFieldVolumeData);

	if (DistanceFieldResolutionScale > 0)
	{
		const double StartTime = FPlatformTime::Seconds();

		FEmbreeScene EmbreeScene;
		SetupEmbreeScene(MeshName,
			SourceMeshData,
			LODModel,
			MaterialBlendModes,
			bGenerateAsIfTwoSided,
			EmbreeScene);

		check(EmbreeScene.bUseEmbree);

		// Whether to use an Embree Point Query to compute the closest unsigned distance.  Rays will only be traced to determine backfaces visible for sign.
		const bool bUsePointQuery = true;

		TArray<FVector3f> SampleDirections;
		{
			const int32 NumVoxelDistanceSamples = bUsePointQuery ? 49 : 576;
			FRandomStream RandomStream(0);
			GenerateStratifiedUniformHemisphereSamples(NumVoxelDistanceSamples, RandomStream, SampleDirections);
			TArray<FVector3f> OtherHemisphereSamples;
			GenerateStratifiedUniformHemisphereSamples(NumVoxelDistanceSamples, RandomStream, OtherHemisphereSamples);

			for (int32 i = 0; i < OtherHemisphereSamples.Num(); i++)
			{
				FVector3f Sample = OtherHemisphereSamples[i];
				Sample.Z *= -1.0f;
				SampleDirections.Add(Sample);
			}
		}

		static const auto CVar = IConsoleManager::Get().FindTConsoleVariableDataInt(TEXT("r.DistanceFields.MaxPerMeshResolution"));
		const int32 PerMeshMax = CVar->GetValueOnAnyThread();

		// Meshes with explicit artist-specified scale can go higher
		const int32 MaxNumBlocksOneDim = FMath::Min<int32>(FMath::DivideAndRoundNearest(DistanceFieldResolutionScale <= 1 ? PerMeshMax / 2 : PerMeshMax, DistanceField::UniqueDataBrickSize), DistanceField::MaxIndirectionDimension - 1);

		static const auto CVarDensity = IConsoleManager::Get().FindTConsoleVariableDataFloat(TEXT("r.DistanceFields.DefaultVoxelDensity"));
		const float VoxelDensity = CVarDensity->GetValueOnAnyThread();

		const float NumVoxelsPerLocalSpaceUnit = VoxelDensity * DistanceFieldResolutionScale;
		FBox LocalSpaceMeshBounds(Bounds.GetBox());

		// Make sure the mesh bounding box has positive extents to handle planes
		{
			FVector MeshBoundsCenter = LocalSpaceMeshBounds.GetCenter();
			FVector MeshBoundsExtent = FVector::Max(LocalSpaceMeshBounds.GetExtent(), FVector(1.0f, 1.0f, 1.0f));
			LocalSpaceMeshBounds.Min = MeshBoundsCenter - MeshBoundsExtent;
			LocalSpaceMeshBounds.Max = MeshBoundsCenter + MeshBoundsExtent;
		}

		// We sample on voxel corners and use central differencing for gradients, so a box mesh using two-sided materials whose vertices lie on LocalSpaceMeshBounds produces a zero gradient on intersection
		// Expand the mesh bounds by a fraction of a voxel to allow room for a pullback on the hit location for computing the gradient.
		// Only expand for two sided meshes as this adds significant Mesh SDF tracing cost
		if (EmbreeScene.bMostlyTwoSided)
		{
			const FVector DesiredDimensions = FVector(LocalSpaceMeshBounds.GetSize() * FVector(NumVoxelsPerLocalSpaceUnit / (float)DistanceField::UniqueDataBrickSize));
			const FIntVector Mip0IndirectionDimensions = FIntVector(
				FMath::Clamp(FMath::RoundToInt(DesiredDimensions.X), 1, MaxNumBlocksOneDim),
				FMath::Clamp(FMath::RoundToInt(DesiredDimensions.Y), 1, MaxNumBlocksOneDim),
				FMath::Clamp(FMath::RoundToInt(DesiredDimensions.Z), 1, MaxNumBlocksOneDim));

			const float CentralDifferencingExpandInVoxels = .25f;
			const FVector TexelObjectSpaceSize = LocalSpaceMeshBounds.GetSize() / FVector(Mip0IndirectionDimensions * DistanceField::UniqueDataBrickSize - FIntVector(2 * CentralDifferencingExpandInVoxels));
			LocalSpaceMeshBounds = LocalSpaceMeshBounds.ExpandBy(TexelObjectSpaceSize);
		}

		// The tracing shader uses a Volume space that is normalized by the maximum extent, to keep Volume space within [-1, 1], we must match that behavior when encoding
		const float LocalToVolumeScale = 1.0f / LocalSpaceMeshBounds.GetExtent().GetMax();

		const FVector DesiredDimensions = FVector(LocalSpaceMeshBounds.GetSize() * FVector(NumVoxelsPerLocalSpaceUnit / (float)DistanceField::UniqueDataBrickSize));
		const FIntVector Mip0IndirectionDimensions = FIntVector(
			FMath::Clamp(FMath::RoundToInt(DesiredDimensions.X), 1, MaxNumBlocksOneDim),
			FMath::Clamp(FMath::RoundToInt(DesiredDimensions.Y), 1, MaxNumBlocksOneDim),
			FMath::Clamp(FMath::RoundToInt(DesiredDimensions.Z), 1, MaxNumBlocksOneDim));

		TArray<uint8> StreamableMipData;

		for (int32 MipIndex = 0; MipIndex < DistanceField::NumMips; MipIndex++)
		{
			const FIntVector IndirectionDimensions = FIntVector(
				FMath::DivideAndRoundUp(Mip0IndirectionDimensions.X, 1 << MipIndex),
				FMath::DivideAndRoundUp(Mip0IndirectionDimensions.Y, 1 << MipIndex),
				FMath::DivideAndRoundUp(Mip0IndirectionDimensions.Z, 1 << MipIndex));

			// Expand to guarantee one voxel border for gradient reconstruction using bilinear filtering
			const FVector TexelObjectSpaceSize = LocalSpaceMeshBounds.GetSize() / FVector(IndirectionDimensions * DistanceField::UniqueDataBrickSize - FIntVector(2 * DistanceField::MeshDistanceFieldObjectBorder));
			const FBox DistanceFieldVolumeBounds = LocalSpaceMeshBounds.ExpandBy(TexelObjectSpaceSize);

			const FVector IndirectionVoxelSize = DistanceFieldVolumeBounds.GetSize() / FVector(IndirectionDimensions);
			const float IndirectionVoxelRadius = IndirectionVoxelSize.Size();

			const FVector VolumeSpaceDistanceFieldVoxelSize = IndirectionVoxelSize * LocalToVolumeScale / FVector(DistanceField::UniqueDataBrickSize);
			const float MaxDistanceForEncoding = VolumeSpaceDistanceFieldVoxelSize.Size() * DistanceField::BandSizeInVoxels;
			const float LocalSpaceTraceDistance = MaxDistanceForEncoding / LocalToVolumeScale;
			const FVector2D DistanceFieldToVolumeScaleBias(2.0f * MaxDistanceForEncoding, -MaxDistanceForEncoding);

			TArray<FSparseMeshDistanceFieldAsyncTask> AsyncTasks;
			AsyncTasks.Reserve(IndirectionDimensions.X * IndirectionDimensions.Y * IndirectionDimensions.Z / 8);

			for (int32 ZIndex = 0; ZIndex < IndirectionDimensions.Z; ZIndex++)
			{
				for (int32 YIndex = 0; YIndex < IndirectionDimensions.Y; YIndex++)
				{
					for (int32 XIndex = 0; XIndex < IndirectionDimensions.X; XIndex++)
					{
						AsyncTasks.Emplace(
							EmbreeScene,
							&SampleDirections,
							LocalSpaceTraceDistance,
							DistanceFieldVolumeBounds,
							LocalToVolumeScale,
							DistanceFieldToVolumeScaleBias,
							FIntVector(XIndex, YIndex, ZIndex),
							IndirectionDimensions,
							bUsePointQuery);
					}
				}
			}

			static bool bMultiThreaded = true;

			if (bMultiThreaded)
			{
				EParallelForFlags Flags = EParallelForFlags::BackgroundPriority | EParallelForFlags::Unbalanced;

				ParallelForTemplate(
					TEXT("GenerateSignedDistanceFieldVolumeData.PF"),
					AsyncTasks.Num(), 1, [&AsyncTasks](int32 TaskIndex)
					{
						AsyncTasks[TaskIndex].DoWork();
					}, Flags);
			}
			else
			{
				for (FSparseMeshDistanceFieldAsyncTask& AsyncTask : AsyncTasks)
				{
					AsyncTask.DoWork();
				}
			}

			FSparseDistanceFieldMip& OutMip = OutData.Mips[MipIndex];
			TArray<uint32> IndirectionTable;
			IndirectionTable.Empty(IndirectionDimensions.X * IndirectionDimensions.Y * IndirectionDimensions.Z);
			IndirectionTable.AddUninitialized(IndirectionDimensions.X * IndirectionDimensions.Y * IndirectionDimensions.Z);

			for (int32 i = 0; i < IndirectionTable.Num(); i++)
			{
				IndirectionTable[i] = DistanceField::InvalidBrickIndex;
			}

			TArray<FSparseMeshDistanceFieldAsyncTask*> ValidBricks;
			ValidBricks.Empty(AsyncTasks.Num());

			for (int32 TaskIndex = 0; TaskIndex < AsyncTasks.Num(); TaskIndex++)
			{
				if (AsyncTasks[TaskIndex].BrickMinDistance < MAX_uint8 && AsyncTasks[TaskIndex].BrickMaxDistance > MIN_uint8)
				{
					ValidBricks.Add(&AsyncTasks[TaskIndex]);
				}
			}

			const uint32 NumBricks = ValidBricks.Num();

			const uint32 BrickSizeBytes = DistanceField::BrickSize * DistanceField::BrickSize * DistanceField::BrickSize * GPixelFormats[DistanceField::DistanceFieldFormat].BlockBytes;

			TArray<uint8> DistanceFieldBrickData;
			DistanceFieldBrickData.Empty(BrickSizeBytes * NumBricks);
			DistanceFieldBrickData.AddUninitialized(BrickSizeBytes * NumBricks);

			for (int32 BrickIndex = 0; BrickIndex < ValidBricks.Num(); BrickIndex++)
			{
				const FSparseMeshDistanceFieldAsyncTask& Brick = *ValidBricks[BrickIndex];
				const int32 IndirectionIndex = ComputeLinearVoxelIndex(Brick.BrickCoordinate, IndirectionDimensions);
				IndirectionTable[IndirectionIndex] = BrickIndex;

				check(BrickSizeBytes == Brick.DistanceFieldVolume.Num() * Brick.DistanceFieldVolume.GetTypeSize());
				FPlatformMemory::Memcpy(&DistanceFieldBrickData[BrickIndex * BrickSizeBytes], Brick.DistanceFieldVolume.GetData(), Brick.DistanceFieldVolume.Num() * Brick.DistanceFieldVolume.GetTypeSize());
			}

			const int32 IndirectionTableBytes = IndirectionTable.Num() * IndirectionTable.GetTypeSize();
			const int32 MipDataBytes = IndirectionTableBytes + DistanceFieldBrickData.Num();

			if (MipIndex == DistanceField::NumMips - 1)
			{
				OutData.AlwaysLoadedMip.Empty(MipDataBytes);
				OutData.AlwaysLoadedMip.AddUninitialized(MipDataBytes);

				FPlatformMemory::Memcpy(&OutData.AlwaysLoadedMip[0], IndirectionTable.GetData(), IndirectionTableBytes);

				if (DistanceFieldBrickData.Num() > 0)
				{
					FPlatformMemory::Memcpy(&OutData.AlwaysLoadedMip[IndirectionTableBytes], DistanceFieldBrickData.GetData(), DistanceFieldBrickData.Num());
				}
			}
			else
			{
				OutMip.BulkOffset = StreamableMipData.Num();
				StreamableMipData.AddUninitialized(MipDataBytes);
				OutMip.BulkSize = StreamableMipData.Num() - OutMip.BulkOffset;
				checkf(OutMip.BulkSize > 0, TEXT("BulkSize was 0 for %s with %ux%ux%u indirection"), *MeshName, IndirectionDimensions.X, IndirectionDimensions.Y, IndirectionDimensions.Z);

				FPlatformMemory::Memcpy(&StreamableMipData[OutMip.BulkOffset], IndirectionTable.GetData(), IndirectionTableBytes);

				if (DistanceFieldBrickData.Num() > 0)
				{
					FPlatformMemory::Memcpy(&StreamableMipData[OutMip.BulkOffset + IndirectionTableBytes], DistanceFieldBrickData.GetData(), DistanceFieldBrickData.Num());
				}
			}

			OutMip.IndirectionDimensions = IndirectionDimensions;
			OutMip.DistanceFieldToVolumeScaleBias = DistanceFieldToVolumeScaleBias;
			OutMip.NumDistanceFieldBricks = NumBricks;

			// Account for the border voxels we added
			const FVector VirtualUVMin = FVector(DistanceField::MeshDistanceFieldObjectBorder) / FVector(IndirectionDimensions * DistanceField::UniqueDataBrickSize);
			const FVector VirtualUVSize = FVector(IndirectionDimensions * DistanceField::UniqueDataBrickSize - FIntVector(2 * DistanceField::MeshDistanceFieldObjectBorder)) / FVector(IndirectionDimensions * DistanceField::UniqueDataBrickSize);

			const FVector VolumePositionExtent = LocalSpaceMeshBounds.GetExtent() * LocalToVolumeScale;

			// [-VolumePositionExtent, VolumePositionExtent] -> [VirtualUVMin, VirtualUVMin + VirtualUVSize]
			OutMip.VolumeToVirtualUVScale = VirtualUVSize / (2 * VolumePositionExtent);
			OutMip.VolumeToVirtualUVAdd = VolumePositionExtent * OutMip.VolumeToVirtualUVScale + VirtualUVMin;
		}

		DeleteEmbreeScene(EmbreeScene);

		OutData.bMostlyTwoSided = EmbreeScene.bMostlyTwoSided;
		OutData.LocalSpaceMeshBounds = LocalSpaceMeshBounds;

		OutData.StreamableMips.Lock(LOCK_READ_WRITE);
		uint8* Ptr = (uint8*)OutData.StreamableMips.Realloc(StreamableMipData.Num());
		FMemory::Memcpy(Ptr, StreamableMipData.GetData(), StreamableMipData.Num());
		OutData.StreamableMips.Unlock();
		OutData.StreamableMips.SetBulkDataFlags(BULKDATA_Force_NOT_InlinePayload);

		const float BuildTime = (float)(FPlatformTime::Seconds() - StartTime);

		if (BuildTime > 1.0f)
		{
			UE_LOG(LogTemp, Log, TEXT("Finished distance field build in %.1fs - %ux%ux%u sparse distance field, %.1fMb total, %.1fMb always loaded, %u%% occupied, %u triangles, %s"),
				BuildTime,
				Mip0IndirectionDimensions.X * DistanceField::UniqueDataBrickSize,
				Mip0IndirectionDimensions.Y * DistanceField::UniqueDataBrickSize,
				Mip0IndirectionDimensions.Z * DistanceField::UniqueDataBrickSize,
				(OutData.GetResourceSizeBytes() + OutData.StreamableMips.GetBulkDataSize()) / 1024.0f / 1024.0f,
				(OutData.AlwaysLoadedMip.GetAllocatedSize()) / 1024.0f / 1024.0f,
				FMath::RoundToInt(100.0f * OutData.Mips[0].NumDistanceFieldBricks / (float)(Mip0IndirectionDimensions.X * Mip0IndirectionDimensions.Y * Mip0IndirectionDimensions.Z)),
				EmbreeScene.NumIndices / 3,
				*MeshName);
		}
	}
}

bool USignedDistanceFieldUtilities::GenerateSDF(UStaticMesh* StaticMesh)
{
	if (!StaticMesh->IsValidLowLevel())
		return false;

	const TArray<FStaticMaterial>& StaticMaterials = StaticMesh->GetStaticMaterials();

	TArray<USignedDistanceFieldUtilities::FSignedDistanceFieldBuildMaterialData> BuildMaterialData;
	BuildMaterialData.SetNum(StaticMaterials.Num());

 	FMeshSectionInfoMap& SectionInfoMap = StaticMesh->GetSectionInfoMap();
	const uint32 LODIndex = 0;

	for (int32 SectionIndex = 0; SectionIndex < SectionInfoMap.GetSectionNumber(LODIndex); SectionIndex++)
	{
		const FMeshSectionInfo& Section = SectionInfoMap.Get(LODIndex, SectionIndex);

		if (!BuildMaterialData.IsValidIndex(Section.MaterialIndex))
		{
			continue;
		}

		USignedDistanceFieldUtilities::FSignedDistanceFieldBuildMaterialData& MaterialData = BuildMaterialData[Section.MaterialIndex];
		MaterialData.bAffectDistanceFieldLighting = Section.bAffectDistanceFieldLighting;

		UMaterialInterface* MaterialInterface = StaticMaterials[Section.MaterialIndex].MaterialInterface;
		if (MaterialInterface)
		{
			MaterialData.BlendMode = MaterialInterface->GetBlendMode();
			MaterialData.bTwoSided = MaterialInterface->IsTwoSided();
		}
	}

	//FString DistanceFieldKey = BuildDistanceFieldDerivedDataKey(InStaticMeshDerivedDataKey);

	//for (int32 MaterialIndex = 0; MaterialIndex < Mesh->GetStaticMaterials().Num(); MaterialIndex++)
	//{
	//	DistanceFieldKey += FString::Printf(TEXT("_M%u_%u_%u"),
	//		(uint32)BuildMaterialData[MaterialIndex].BlendMode,
	//		BuildMaterialData[MaterialIndex].bTwoSided ? 1 : 0,
	//		BuildMaterialData[MaterialIndex].bAffectDistanceFieldLighting ? 1 : 0);
	//}

	FString MeshName = StaticMesh->GetName();

	const FMeshBuildSettings& BuildSettings = StaticMesh->GetSourceModel(0).BuildSettings;
	UStaticMesh* GenerateSource = BuildSettings.DistanceFieldReplacementMesh ? ToRawPtr(BuildSettings.DistanceFieldReplacementMesh) : StaticMesh;
	float DistanceFieldResolutionScale = BuildSettings.DistanceFieldResolutionScale;
	bool bGenerateDistanceFieldAsIfTwoSided = BuildSettings.bGenerateDistanceFieldAsIfTwoSided;

	FDistanceFieldVolumeData* GeneratedVolumeData = new FDistanceFieldVolumeData();
	FSourceMeshDataForDerivedDataTask SourceMeshData{};
	if (GenerateSource->GetRenderData())
	{
		const FStaticMeshLODResources& LODModel = GenerateSource->GetRenderData()->LODResources[0];

		//USignedDistanceFieldUtilities* MyClass = NewObject<USignedDistanceFieldUtilities>();
		/*MyClass->*/GenerateSignedDistanceFieldVolumeData(
			MeshName,
			SourceMeshData,
			LODModel,
			MoveTemp(BuildMaterialData),
			GenerateSource->GetRenderData()->Bounds,
			DistanceFieldResolutionScale,
			bGenerateDistanceFieldAsIfTwoSided,
			*GeneratedVolumeData
		);


		// Editor 'force delete' can null any UObject pointers which are seen by reference collecting (eg FProperty or serialized)
		//if (Task->StaticMesh)
		{
			FObjectCacheContextScope ObjectCacheScope;

			check(!StaticMesh->IsCompiling());

			GeneratedVolumeData->bAsyncBuilding = false;

			FStaticMeshRenderData* RenderData = StaticMesh->GetRenderData();
			FDistanceFieldVolumeData* OldVolumeData = RenderData->LODResources[0].DistanceFieldData;

			// Assign the new volume data, this is safe because the render thread makes a copy of the pointer at scene proxy creation time.
			RenderData->LODResources[0].DistanceFieldData = GeneratedVolumeData;

			// Renderstates are not initialized between UStaticMesh::PreEditChange() and UStaticMesh::PostEditChange()
			if (RenderData->IsInitialized())
			{
				for (UStaticMeshComponent* Component : ObjectCacheScope.GetContext().GetStaticMeshComponents(StaticMesh))
				{
					if (Component->IsRegistered() && Component->IsRenderStateCreated())
					{
						Component->MarkRenderStateDirty();
					}
				}
			}

			if (OldVolumeData)
			{
				// Rendering thread may still be referencing the old one, use the deferred cleanup interface to delete it next frame when it is safe
				BeginCleanup(OldVolumeData);
			}

			// Need also to update platform render data if it's being cached
			FStaticMeshRenderData* PlatformRenderData = RenderData->NextCachedRenderData.Get();
			while (PlatformRenderData)
			{
				if (PlatformRenderData->LODResources[0].DistanceFieldData)
				{
					*PlatformRenderData->LODResources[0].DistanceFieldData = *GeneratedVolumeData;
					// The old bulk data assignment operator doesn't copy over flags
					PlatformRenderData->LODResources[0].DistanceFieldData->StreamableMips.ResetBulkDataFlags(GeneratedVolumeData->StreamableMips.GetBulkDataFlags());
				}
				PlatformRenderData = PlatformRenderData->NextCachedRenderData.Get();
			}

			//{
			//	TArray<uint8> DerivedData;
			//	// Save built distance field volume to DDC
			//	FMemoryWriter Ar(DerivedData, /*bIsPersistent=*/ true);
			//	StaticMesh->GetRenderData()->LODResources[0].DistanceFieldData->Serialize(Ar, Task->StaticMesh);
			//	GetDerivedDataCacheRef().Put(*Task->DDCKey, DerivedData, Task->StaticMesh->GetPathName());
			//	COOK_STAT(Timer.AddMiss(DerivedData.Num()));
			//}

			//BeginCacheMeshCardRepresentation(
			//	Task->TargetPlatform,
			//	Task->StaticMesh,
			//	Task->StaticMesh->GetPlatformStaticMeshRenderData(Task->StaticMesh, Task->TargetPlatform),
			//	Task->DDCKey,
			//	&Task->SourceMeshData);
		}

		return true;
	}

	return false;
}

#else
//
//void FMeshUtilities::GenerateSignedDistanceFieldVolumeData(
//	FString MeshName,
//	const FSourceMeshDataForDerivedDataTask& SourceMeshData,
//	const FStaticMeshLODResources& LODModel,
//	class FQueuedThreadPool& ThreadPool,
//	const TArray<FSignedDistanceFieldBuildMaterialData>& MaterialBlendModes,
//	const FBoxSphereBounds& Bounds,
//	float DistanceFieldResolutionScale,
//	bool bGenerateAsIfTwoSided,
//	FDistanceFieldVolumeData& OutData)
//{
//	if (DistanceFieldResolutionScale > 0)
//	{
//		UE_LOG(LogTemp, Warning, TEXT("Couldn't generate distance field for mesh, platform is missing Embree support."));
//	}
//}

#endif // PLATFORM_ENABLE_VECTORINTRINSICS


