#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Conservative.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Fp.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_ReactionMap.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Reactive.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Return.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Statistics.cpp
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" BeadPack_Strips.cpp
	rm *.bak
	exit 0
fi
if [[ "$#" -eq 2 ]]
then
	sed -i'.bak' "s/.*using namespace beadpack::model_.*/  using namespace beadpack::model_$1;/" "BeadPack_$2.cpp"
	rm *.bak
	exit 0
fi
echo "Bad parameters."
exit 1
