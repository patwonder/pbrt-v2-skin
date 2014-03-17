// MultipoleProfileCalculator.cpp : 定义 DLL 应用程序的导出函数。
//

#include "stdafx.h"
#include "MultipoleProfileCalculator.h"


// 这是导出变量的一个示例
MULTIPOLEPROFILECALCULATOR_API int nMultipoleProfileCalculator=0;

// 这是导出函数的一个示例。
MULTIPOLEPROFILECALCULATOR_API int fnMultipoleProfileCalculator(void)
{
	return 42;
}

// 这是已导出类的构造函数。
// 有关类定义的信息，请参阅 MultipoleProfileCalculator.h
CMultipoleProfileCalculator::CMultipoleProfileCalculator()
{
	return;
}
