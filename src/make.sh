#!/bin/bash
if [[ "$#" -eq 1 ]]
then
	./set_model.sh $1
	make all
	mv ../bin/BeadPack_Conservative "../bin/BeadPack_Conservative_$1"
	mv ../bin/BeadPack_Fp "../bin/BeadPack_Fp_$1"
	mv ../bin/BeadPack_ReactionMap "../bin/BeadPack_ReactionMap_$1"
	mv ../bin/BeadPack_Reactive "../bin/BeadPack_Reactive_$1"
	mv ../bin/BeadPack_Return "../bin/BeadPack_Return_$1"
	mv ../bin/BeadPack_Statistics "../bin/BeadPack_Statistics_$1"
	mv ../bin/BeadPack_Strips "../bin/BeadPack_Strips_$1"
  mv ../bin/BeadPack_Poincare "../bin/BeadPack_Poincare_$1"
	exit 0
fi
if [[ "$#" -eq 2 ]]
then
	./set_model.sh $1 $2
	make "BeadPack_$2"
	mv "../bin/BeadPack_$2" "../bin/BeadPack_$2_$1"
	exit 0
fi
echo "Bad parameters."
exit 1
