// Copyright Epic Games, Inc. All Rights Reserved.

#pragma once

#include "CoreMinimal.h"
#include "GameFramework/GameModeBase.h"
#include "DistanceFieldTestGameModeBase.generated.h"

/**
 * 
 */
UCLASS()
class DISTANCEFIELDTEST_API ADistanceFieldTestGameModeBase : public AGameModeBase
{
	GENERATED_BODY()

public:

    ADistanceFieldTestGameModeBase();

protected:

    virtual void BeginPlay() override;

private:

    //UMG蓝图的实例对象，用于显示在游戏的Viewport中
    UUserWidget* MyWidgetInstance;
};
