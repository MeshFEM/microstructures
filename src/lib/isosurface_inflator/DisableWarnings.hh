// Disable compiler warnings before including 3rdparty code
#if defined(__clang__)
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wshadow"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wsign-compare"
	#pragma clang diagnostic push
	#pragma clang diagnostic ignored "-Wswitch-default"
#elif (defined(__GNUC__) || defined(__GNUG__)) && !(defined(__clang__) || defined(__INTEL_COMPILER))
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wshadow"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare"
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wswitch-default"
#endif
