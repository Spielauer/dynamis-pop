#include "model.h"

namespace mm {


double	Tab0Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab0Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab0Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab1Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab1Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab1Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab2Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(1, nIndex);
}

double	Tab2Expr1( XTable *prTable, int nIndex )
{
	return prTable->GetVal(2, nIndex);
}

double	Tab2Expr2( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab3Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab4Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(2, nIndex);
}

double	Tab4Expr1( XTable *prTable, int nIndex )
{
	return prTable->GetVal(1, nIndex);
}

double	Tab4Expr2( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab5Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(3, nIndex);
}

double	Tab5Expr1( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab5Expr2( XTable *prTable, int nIndex )
{
	return prTable->GetVal(1, nIndex);
}

double	Tab5Expr3( XTable *prTable, int nIndex )
{
	return prTable->GetVal(2, nIndex);
}

double	Tab5Expr4( XTable *prTable, int nIndex )
{
	return prTable->GetVal(4, nIndex);
}

double	Tab6Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab7Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab7Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab7Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab8Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab8Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab8Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab9Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab10Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab11Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab12Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab13Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab13Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab13Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab14Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(2, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab14Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab15Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(5, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab15Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(7, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab15Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab15Expr3( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(6, nIndex), prTable->GetVal(2, nIndex) );
}

double	Tab15Expr4( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(8, nIndex), prTable->GetVal(2, nIndex) );
}

double	Tab15Expr5( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(4, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab16Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(5, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab16Expr1( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(7, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab16Expr2( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(3, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab16Expr3( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(6, nIndex), prTable->GetVal(2, nIndex) );
}

double	Tab16Expr4( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(8, nIndex), prTable->GetVal(2, nIndex) );
}

double	Tab16Expr5( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(4, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab17Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab18Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab19Expr0( XTable *prTable, int nIndex )
{
	return prTable->GetVal(0, nIndex);
}

double	Tab20Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(0, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab21Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(0, nIndex), prTable->GetVal(1, nIndex) );
}

double	Tab22Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab23Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

double	Tab24Expr0( XTable *prTable, int nIndex )
{
	return Division( prTable->GetVal(1, nIndex), prTable->GetVal(0, nIndex) );
}

#pragma optimize( "", off )
void InitTables()
{
	AddTable( 0,
		Tab0Expr0,
		Tab0Expr1,
		Tab0Expr2 );
	AddTable( 1,
		Tab1Expr0,
		Tab1Expr1,
		Tab1Expr2 );
	AddTable( 2,
		Tab2Expr0,
		Tab2Expr1,
		Tab2Expr2 );
	AddTable( 3,
		Tab3Expr0 );
	AddTable( 4,
		Tab4Expr0,
		Tab4Expr1,
		Tab4Expr2 );
	AddTable( 5,
		Tab5Expr0,
		Tab5Expr1,
		Tab5Expr2,
		Tab5Expr3,
		Tab5Expr4 );
	AddTable( 6,
		Tab6Expr0 );
	AddTable( 7,
		Tab7Expr0,
		Tab7Expr1,
		Tab7Expr2 );
	AddTable( 8,
		Tab8Expr0,
		Tab8Expr1,
		Tab8Expr2 );
	AddTable( 9,
		Tab9Expr0 );
	AddTable( 10,
		Tab10Expr0 );
	AddTable( 11,
		Tab11Expr0 );
	AddTable( 12,
		Tab12Expr0 );
	AddTable( 13,
		Tab13Expr0,
		Tab13Expr1,
		Tab13Expr2 );
	AddTable( 14,
		Tab14Expr0,
		Tab14Expr1 );
	AddTable( 15,
		Tab15Expr0,
		Tab15Expr1,
		Tab15Expr2,
		Tab15Expr3,
		Tab15Expr4,
		Tab15Expr5 );
	AddTable( 16,
		Tab16Expr0,
		Tab16Expr1,
		Tab16Expr2,
		Tab16Expr3,
		Tab16Expr4,
		Tab16Expr5 );
	AddTable( 17,
		Tab17Expr0 );
	AddTable( 18,
		Tab18Expr0 );
	AddTable( 19,
		Tab19Expr0 );
	AddTable( 20,
		Tab20Expr0 );
	AddTable( 21,
		Tab21Expr0 );
	AddTable( 22,
		Tab22Expr0 );
	AddTable( 23,
		Tab23Expr0 );
	AddTable( 24,
		Tab24Expr0 );
}

#pragma optimize( "", on )


int nTesTPaths = CheckPathsForLibraries_Modgen12();

} // namespace 
