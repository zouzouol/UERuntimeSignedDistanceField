// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/HUD.h"
#include "MyHUD.generated.h"

/**
 * 
 */
UCLASS()
class DISTANCEFIELDTEST_API AMyHUD : public AHUD
{
	GENERATED_BODY()
public:

	UPROPERTY(BlueprintReadWrite, EditAnywhere)
	TSubclassOf<UUserWidget> MainWidgetClass;
	UPROPERTY(BlueprintReadWrite, EditAnywhere)
	TSubclassOf<UUserWidget> StartWidgetClass;
	virtual void BeginPlay() override;

	virtual void DrawHUD() override;
	bool bDrawDebugText;
	TArray<FString> DebugTexts;
};
