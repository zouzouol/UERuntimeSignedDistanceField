// Fill out your copyright notice in the Description page of Project Settings.


#include "MyHUD.h"

#include "Blueprint/WidgetBlueprintLibrary.h"
#include "Blueprint/UserWidget.h"
#include "Engine/Canvas.h"

void AMyHUD::BeginPlay()
{
	if (StartWidgetClass)
	{
		if (UUserWidget* MainWidget = CreateWidget<UUserWidget>(GetWorld()->GetFirstPlayerController(), StartWidgetClass))
		{
			MainWidget->AddToViewport();
			//if (ANewProjectPlayerController* PlayerController = Cast<ANewProjectPlayerController>(GetWorld()->GetFirstPlayerController()))
			//{
			//	PlayerController->StartWidget = MainWidget;
			//	return;
			//}
			//UE_LOG(LogTemp, Error, TEXT("����StartWidgetʧ��,�޷��ҵ�������"));
		}
	}
	else
	{
		if (MainWidgetClass)
		{
			if (UUserWidget* MainWidget = CreateWidget<UUserWidget>(GetWorld()->GetFirstPlayerController(), MainWidgetClass))
			{
				MainWidget->AddToViewport();
				//if (ANewProjectPlayerController* PlayerController = Cast<ANewProjectPlayerController>(GetWorld()->GetFirstPlayerController()))
				//{
				//	PlayerController->MainWidget = MainWidget;
				//	PlayerController->InitAfterHUD();
				//	return;
				//}
				//UE_LOG(LogTemp, Error, TEXT("����UMG_Mainʧ��,�޷��ҵ�������"));
				return;
			}
		}
		//UE_LOG(LogTemp, Error, TEXT("����UMG_Mainʧ��,�޷��ҵ���ͼ·��"));
	}
}

void AMyHUD::DrawHUD()
{
	Super::DrawHUD();


	if (bDrawDebugText)
	{
		for (int32 i = 0; i < DebugTexts.Num(); ++i)
		{
			DrawText(FString::Printf(TEXT("%s"), *DebugTexts[i]), FColor::White, Canvas.Get()->SizeX - 128,
				Canvas.Get()->SizeY - (i + 1) * 32, 0, 1.5f);
		}
	}
}
