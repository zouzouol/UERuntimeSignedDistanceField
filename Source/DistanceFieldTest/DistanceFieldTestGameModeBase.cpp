// Copyright Epic Games, Inc. All Rights Reserved.


#include "DistanceFieldTestGameModeBase.h"
#include "Blueprint/UserWidget.h"

ADistanceFieldTestGameModeBase::ADistanceFieldTestGameModeBase()
{
}

void ADistanceFieldTestGameModeBase::BeginPlay()
{
    //���Widget�����Ƿ���ڣ�����������Ƴ�����
    if (MyWidgetInstance)
    {
        MyWidgetInstance->RemoveFromViewport();
        MyWidgetInstance = nullptr;
    }

    ///Script/UMGEditor.WidgetBlueprint'/Game/WBP_NewWidget2.WBP_NewWidget2'
    ///Script/UMGEditor.WidgetBlueprint'/Game/BPW_MainWin.BPW_MainWin'
    ///Script/UMGEditor.WidgetBlueprint'/Game/WBP_Main.WBP_Main'
    // 
    
    ////�����Զ���UMG��class��ͨ�����class����Widget���󣬲���ʾ�ڽ����С�
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
