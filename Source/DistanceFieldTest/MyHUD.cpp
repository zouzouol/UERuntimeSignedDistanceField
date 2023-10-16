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
			//UE_LOG(LogTemp, Error, TEXT("加载StartWidget失败,无法找到控制器"));
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
				//UE_LOG(LogTemp, Error, TEXT("加载UMG_Main失败,无法找到控制器"));
				return;
			}
		}
		//UE_LOG(LogTemp, Error, TEXT("加载UMG_Main失败,无法找到蓝图路径"));
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
