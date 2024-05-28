      module mulliken_mod
!
!     This is a module file to support Mulliken population analysis.
!
!
!     H. P. Hratchian
!     Department of Chemistry & Chemical Biology
!     Center for Chemical Computation and Theory
!     University of California, Merced
!     hhratchian@ucmerced.edu
!
!
!
!     USE Connections
!
      use mqc_general
      use mqc_gaussian
      use mqc_algebra2
      use iso_fortran_env
!
!     Variable Declarations
!
      implicit none
      integer,parameter::IOut=6
!
!
      CONTAINS

!
!PROCEDURE mulliken_analysis_ao
      subroutine mulliken_analysis_ao(PMatrix,SMatrix,mullikenMatrixAO,  &
        mullikenGrossPopAO,basisFunctionTypeList)
!
!     This subroutine carries out Mulliken analysis. The atomic orbital basis is
!     assumed. Input arguments are the density matrix (<PMatrix>) and the AO
!     overlap matrix (<SMatrix>). Optional output arguments are the full
!     Mulliken matrix in the AO basis (<mullikenMatrixAO>) and gross populations
!     (<mullikenGrossPopAO>).
!
!
!     - H. P. Hratchian, 2024.
!
!
!     Variable Declarations
!
      type(MQC_Variable),intent(in)::PMatrix,SMatrix
      type(MQC_Variable),optional::mullikenMatrixAO,mullikenGrossPopAO,  &
        basisFunctionTypeList
      integer(kind=int64)::i,tmpInteger,angularMomentum
      real(kind=real64)::tmpReal
      real(kind=real64),dimension(10)::myMullikenGrossPopAngMomentum
      type(MQC_Variable)::myMullikenMatrixAO,myMullikenGrossPopAO
!
!     Form the Mulliken matrix and compress it down to gross populations.
!
      myMullikenMatrixAO = PMatrix*SMatrix
      myMullikenGrossPopAO = MQC_Variable_Array_Sum(myMullikenMatrixAO,1)
!
!     If requested, loop through the gross population values and aggregate gross
!     populations by angular momentum.
!
      write(*,*)
      call myMullikenGrossPopAO%print(iOut,header='grossPops')
      write(*,*)
      if(PRESENT(basisFunctionTypeList)) then
        myMullikenGrossPopAngMomentum = 0
        do i = 1, SIZE(myMullikenGrossPopAO)
          tmpReal = myMullikenGrossPopAO%getVal([i])
          tmpInteger = basisFunctionTypeList%getVal([i])
          call basisType2FunctionInfo(tmpInteger,angularMomentum=angularMomentum)
          write(*,*)' MO i: ',i
          write(*,*)'        population val  : ',tmpReal
          write(*,*)'        basis type value: ',tmpInteger
          write(*,*)'        angular momentum: ',angularMomentum
          write(*,*)
          myMullikenGrossPopAngMomentum(angularMomentum+1) =  &
            myMullikenGrossPopAngMomentum(angularMomentum+1) + tmpReal
        endDo
        call mqc_print(myMullikenGrossPopAngMomentum,iOut,header='pops by angular mom')
      endIf
!
!     Fill the possible return arguments that the calling program unit
!     requested.
!
      if(PRESENT(mullikenMatrixAO)) mullikenMatrixAO = myMullikenMatrixAO
      if(PRESENT(mullikenGrossPopAO)) mullikenGrossPopAO = myMullikenGrossPopAO
!
      return
      end subroutine mulliken_analysis_ao


!
!
!
      end module mulliken_mod
