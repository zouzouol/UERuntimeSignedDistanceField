// Fill out your copyright notice in the Description page of Project Settings.

#pragma once

#include "CoreMinimal.h"
#include "Blueprint/UserWidget.h"
#include "MainWidget.generated.h"

/**
 * 
 */
UCLASS()
class DISTANCEFIELDTEST_API UMainWidget : public UUserWidget
{
	GENERATED_BODY()
public:
	UMainWidget(const FObjectInitializer& ObjectInitializer);

protected:
	virtual void NativeConstruct() override;

	UFUNCTION()
	void ButtonClicked();

	UFUNCTION()
	void ButtonSDFClicked();
};
