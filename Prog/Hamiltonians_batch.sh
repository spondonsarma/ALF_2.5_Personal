#!/bin/sh
echo "interface" > "Hamiltonians_interface.h"
rm -f "Hamiltonians_case.h"

while read -r ham
do
  echo "${ham}"
  echo "  module subroutine Ham_Set_${ham}()" >> "Hamiltonians_interface.h"
  echo "  end subroutine Ham_Set_${ham}" >> "Hamiltonians_interface.h"

  echo "Case ('${ham}')" >> "Hamiltonians_case.h"
  echo "  call Ham_Set_${ham}()" >> "Hamiltonians_case.h"
done < Hamiltonians.list

echo "end interface" >> "Hamiltonians_interface.h"
