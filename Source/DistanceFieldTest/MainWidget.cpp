// Fill out your copyright notice in the Description page of Project Settings.


#include "MainWidget.h"

#include "SignedDistanceFieldUtilities.h"

#include "Engine/StaticMesh.h"
#include "Engine/StaticMeshActor.h"
#include "Kismet/GameplayStatics.h"
#include "Components/Button.h"
#include "DatasmithRuntime.h"
#include "DatasmithRuntimeBlueprintLibrary.h"
//#include "DistanceFieldAtlas.h"
//#include "ObjectCacheContext.h"

#include "StaticMeshAttributes.h"
#include "Engine/StaticMeshSourceData.h"

#include "DatasmithAssetUserData.h"


#define MATERIAL_PATH_WHITE TEXT("Material'/Game/Program/Material/M_baise_001.M_baise_001'")
#define MATERIAL_PATH_GLASS TEXT("Material'/Game/Program/Material/M_boli_001.M_boli_001'")
#define CUSTOM_ELEMENT TEXT("CustomElement")
#define ELECTROMECHANICAL TEXT("机电")
#define PHASE_CREATION TEXT("Element*创建的阶段")
#define SAVE_GAME_NAME TEXT("NSSTORAGE")
#define TYPE_GLASS TEXT("Type*玻璃材质")
#define VIRTUAL_EL TEXT("Element*虚拟EL")


#define WALLS TEXT("基本墙")

ADatasmithRuntimeActor* Anchor;

UMainWidget::UMainWidget(const FObjectInitializer& ObjectInitializer) : Super(ObjectInitializer)
{
}

void UMainWidget::NativeConstruct()
{
	Super::NativeConstruct();

	FTransform UserTransformPtr;
	UserTransformPtr.SetLocation(FVector(200.0f, 200.0f, 0.0f));
	FActorSpawnParameters Param;
	Param.SpawnCollisionHandlingOverride = ESpawnActorCollisionHandlingMethod::AlwaysSpawn;
	Anchor = Cast<ADatasmithRuntimeActor>(GetWorld()->SpawnActor(ADatasmithRuntimeActor::StaticClass(), &UserTransformPtr));

    //根据组件ID查找Button组件，并为其添加Click回调事件
    if (UButton* btn = Cast<UButton>(GetWidgetFromName("Button_63")))
    {
        FScriptDelegate Del;
        Del.BindUFunction(this, "ButtonClicked");
        btn->OnClicked.Add(Del);
    }

    if (UButton* btn = Cast<UButton>(GetWidgetFromName("ButtonSDF")))
    {
        FScriptDelegate Del;
        Del.BindUFunction(this, "ButtonSDFClicked");
        btn->OnClicked.Add(Del);
    }
}

void UMainWidget::ButtonClicked()
{
    UE_LOG(LogTemp, Log, TEXT("ButtonClicked"));
    FString DefaultPath = TEXT("D:\\0RevitFiles\\");
    //UDatasmithRuntimeLibrary::LoadFileFromExplorer(Anchor, DefaultPath);


	FString FilePath = TEXT("D:\\datasmithTest\\tttest.udatasmith");
	UDatasmithRuntimeLibrary::LoadFile(Anchor, FilePath);
}

void UMainWidget::ButtonSDFClicked()
{
	TArray<AActor*> OutMeshActor;
	UGameplayStatics::GetAllActorsOfClass(GetWorld(), AActor::StaticClass(), OutMeshActor);
	UE_LOG(LogTemp, Log, TEXT("start Generate SDF data"));
	GEngine->AddOnScreenDebugMessage(0, 400.0f, FColor::Red, TEXT("start Generate SDF data"));

	TSet<UStaticMesh*> StaticMeshSet;
	int32 Index = 0;
	for (AActor* OutActor : OutMeshActor)
	{
		int tagNum = OutActor->Tags.Num();
		if (tagNum > 0)
		{

			UDatasmithAssetUserData* DatasmithAssetUserData = OutActor->GetRootComponent()->GetAssetUserData<UDatasmithAssetUserData>();
			if (DatasmithAssetUserData)
			{

				const FString CreationPhaseKey = PHASE_CREATION;
				const FString CreationPhase = *DatasmithAssetUserData->MetaData.FindRef(*CreationPhaseKey);
				if (CreationPhase == ELECTROMECHANICAL)
				{
					continue;
				}

				//FString familyName = *DatasmithAssetUserData->MetaData.FindRef(TEXT("Element*Family"));
				//if (familyName == TEXT("基本墙") || familyName == TEXT("天花板") || familyName == TEXT("楼板"))
				{
					UStaticMesh* StaticMesh = nullptr;
					if (AStaticMeshActor* StaticMeshActor = Cast<AStaticMeshActor>(OutActor))
					{
						//StaticMeshActor->SetMobility(EComponentMobility::Static);

						auto component = StaticMeshActor->GetStaticMeshComponent();

						component->UnregisterComponent();
						component->RegisterComponent();
						StaticMesh = component->GetStaticMesh();
						if (StaticMesh->IsValidLowLevel())
						{
							StaticMeshSet.Emplace(StaticMesh);
							//break;
						}
					}
				}

				//UE_LOG(LogTemp, Log, TEXT(" familyName555 : %s "), *familyName);
			}


		}
	}

	static int i;

	USignedDistanceFieldUtilities* MyClass = NewObject<USignedDistanceFieldUtilities>();
	for (UStaticMesh* Mesh : StaticMeshSet)
	{
		MyClass->GenerateSDF(Mesh);
		//break;
		i++;
	}

	UE_LOG(LogTemp, Log, TEXT("end Generate SDF data"));

	UE_LOG(LogTemp, Log, TEXT(" Generate SDF data count : %d "), i);


	GEngine->AddOnScreenDebugMessage(0, 400.0f, FColor::Red, TEXT("end Generate SDF data"));
}

