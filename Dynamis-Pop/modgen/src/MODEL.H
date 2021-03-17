#define TIME_TYPE_DOUBLE

//If following file is not found by compiler, check paths for include file and for libraries.
//For information on how to do that, please consult developer's guide.
//Si le fichier suivant n'est pas trouvé par le compilateur, vérifiez les chemins pour les fichiers inclus et les librairies.
//Pour savoir comment modifier ces chemins, veuillez consulter le guide du concepteur de Modgen.
#include "Modgen12.h"

#include "stdafx.h"
#include "common.h"
namespace mm {
	typedef class Calibrator* Calibrator_ptr;
	typedef class Clock* Clock_ptr;
	typedef class Observation* Observation_ptr;
	typedef class Person* Person_ptr;
	typedef class StartpopValues* StartpopValues_ptr;

	inline void initialize_attributes() {}
	inline void enter_simulation() {}
	inline void exit_simulation() {}

	inline void hook_BirthdayEvent() {}
	inline void hook_MidYear() {}
	inline void hook_YearEnd() {}
	inline void hook_YearStart() {}
}
#include "custom_early.h"
#include "actors.h"
#include "tabinit.h"
#include "custom.h"
	
//
// The following code was inserted by patch_modgenXX_outputs after Modgen compilation:
//

// minimum of two doubles
inline double mind(double a, double b)
{
    return (b < a ? b : a);
}

// maximum of two doubles
inline double maxd(double a, double b)
{
    return (a < b ? b : a);
}

// clamped value of a double
inline double clampd(double v, double lo, double hi)
{
    return (v < lo ? lo : v > hi ? hi : v);
}

// minimum of two ints
inline int mini(int a, int b)
{
    return (b < a ? b : a);
}

// maximum of two ints
inline int maxi(int a, int b)
{
    return (a < b ? b : a);
}

// clamped value of an int
inline int clampi(int v, int lo, int hi)
{
    return (v < lo ? lo : v > hi ? hi : v);
}
