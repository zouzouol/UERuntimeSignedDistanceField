// Copyright Epic Games, Inc. All Rights Reserved.


#include "DistanceFieldTestGameModeBase.h"
#include "Blueprint/UserWidget.h"

ADistanceFieldTestGameModeBase::ADistanceFieldTestGameModeBase()
{
}

void ADistanceFieldTestGameModeBase::BeginPlay()
{
    //检测Widget对象是否存在，如果存在则移除掉。
    if (MyWidgetInstance)
    {
        MyWidgetInstance->RemoveFromViewport();
        MyWidgetInstance = nullptr;
    }

    ///Script/UMGEditor.WidgetBlueprint'/Game/WBP_NewWidget2.WBP_NewWidget2'
    ///Script/UMGEditor.WidgetBlueprint'/Game/BPW_MainWin.BPW_MainWin'
    ///Script/UMGEditor.WidgetBlueprint'/Game/WBP_Main.WBP_Main'
    // 
    
    ////加载自定义UMG的class，通过这个class创建Widget对象，并显示在界面中。
    //if (UClass* MyWidgetClass = LoadClass<UUserWidget>(GetWorld()->GetFirstPlayerController(), TEXT("WidgetBlueprint'/Game/BPW_MainWin.BPW_MainWin_C'")))
    //{
    //    if (APlayerController* PC = GetWorld()->GetFirstPlayerController())
    //    {
    //        MyWidgetInstance = CreateWidget<UUserWidget>(PC, MyWidgetClass);
    //        if (MyWidgetInstance)
    //        {
    //            MyWidgetInstance->AddToViewport();
    //        }
    //    }
    //}
}
