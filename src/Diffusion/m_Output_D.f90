module m_output_D
!
! Purpose:
! This module contains the subroutines for printing the flux, k_eff and neutron balance
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties_D
use m_Constants_mod
use m_Results_D
implicit none

contains
    
    subroutine print_flux_diffusion(Prop, Results, N, this)

        type(PropertiesTypeD)   :: Prop
        type(ResultsTypeD)      :: Results
        type(NTypeD)            :: N
        type(Mesh)              :: this

        integer :: j,k, top_right_loc
        real(kind = 8) :: origin(N%D)
        real(kind = 8) :: top_right_coord(N%D)

        top_right_loc = 1

        do j = 1,N%N
            if(N%D == 1) then
                origin = [Prop%Left_B]
                top_right_coord = [Prop%Right_B]
                if (this%Nodes(j,1) == Prop%Right_B) top_right_loc = j
            end if
            if(N%D == 2) then
                origin = [Prop%Left_B, Prop%Bottom_B]
                top_right_coord = [Prop%Right_B, Prop%Top_B]
                if (this%Nodes(j,1) == Prop%Right_B .and. this%Nodes(j,2) == Prop%Top_B) top_right_loc = j
            end if
            if(N%D == 3) then
                origin = [Prop%Left_B, Prop%Bottom_B, Prop%Front_B]
                top_right_coord = [Prop%Right_B, Prop%Top_B, Prop%Back_B]
                if (this%Nodes(j,1) == Prop%Right_B .and. this%Nodes(j,2) == Prop%Top_B .and. this%Nodes(j,3) == Prop%Back_B) top_right_loc = j
            end if
        end do

        print *, "*** Response Results ***"
        print *, " "
        write(*,110) " Position:", [origin]
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Flux(k,1)
        end do
        print *, " "
        print *, "Total ", sum(Results%Flux(:,1))
        print *, " "
        write(*,110) " Position:", this%Nodes((MAXLOC(Results%Flux(1,:))),:)
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Flux(k,(MAXLOC(Results%Flux(1,:))))
        end do
        print *, " "
        print *, "Total ", sum(Results%Flux(:,(MAXLOC(Results%Flux(1,:)))))
        print *, " "
        write(*,110) " Position:" , top_right_coord
        110 format (A,F8.4,F8.4,F8.4)
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Flux(k,top_right_loc)
        end do
        print *, " "
        print *, "Total ", sum(Results%Flux(:,top_right_loc))
        print *, " "

    end subroutine print_flux_diffusion

    subroutine neutron_balance_diffusion(Properties, Results, N)

        type(PropertiesTypeD)   :: Properties
        type(ResultsTypeD)      :: Results
        type(NTypeD)            :: N

        real(kind=8) :: fission_source = 0.0_8, source_total = 0.0_8
        real(kind=8) :: removed = 0.0_8, removed_total = 0.0_8, scattered_in = 0.0_8
        real(kind=8) :: escaped = 0.0_8, escaped_total = 0.0_8, incoming_total = 0.0_8
        real(kind=8) :: alpha = 0.0_8
        integer :: i, j, k

        print *, "Neutron Balance:"
        print *, " "

        if (Properties%LBC == 3) alpha = Properties%Alpha

        do k = 1, N%Group

            fission_source = 0.0_8
            removed = 0.0_8
            scattered_in = 0.0_8
            escaped = 0.0_8

            print *, " Group", k
            print *, " "

            do i = 1, N%Element

                do j = 1, N%Group

                    fission_source = fission_source + (1/Results%k_eff)*Properties%Chi(k)*Properties%Elements(i)%Sigma_f(j)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(j,:)))

                end do

                do j = k,2,-1 ! Upscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(j-1,k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(j-1,:)))

                end do

                do j = k,N%Group-1 ! Downscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(j+1,k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(j+1,:)))

                end do

                removed = removed + Properties%Elements(i)%Sigma_r(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(k,:)))

                escaped = escaped + sum(matmul(((1.0_8/(3.0_8*Properties%Elements(i)%Sigma_t(k)))*Properties%Elements(i)%D_Matrix + Properties%Elements(i)%B_Matrix),Properties%Elements(i)%Flux(k,:)))

            end do

            if (Properties%g == 1) then

                fission_source = fission_source*2.0_8*PI
                scattered_in = scattered_in*2.0_8*PI
                removed = removed*2.0_8*PI
                escaped = escaped*2.0_8*PI

            else if (Properties%g == 2) then

                fission_source = fission_source*4.0_8*PI
                scattered_in = scattered_in*4.0_8*PI
                removed = removed*4.0_8*PI
                escaped = escaped*4.0_8*PI

            end if

            source_total = source_total + fission_source + scattered_in
            print *, "  Source Input      = ", fission_source + scattered_in

            removed_total = removed_total + removed
            print *, "  Removal Rate      = ", removed

            print *, "  Incoming Current  = ", escaped*alpha/(1+alpha)

            if (Properties%LBC == 3) then
                incoming_total = incoming_total + escaped*alpha/(1+alpha)
                escaped = escaped/(1+alpha)
            end if

            escaped_total = escaped_total + escaped
            print *, "  Outgoing Current  = ", escaped

            print *, "  Excess Inflow     = ", fission_source + scattered_in - escaped - removed
            print *, " "

        end do

        if (N%Group > 1)  then

            print *, "Total Source Input      = ", source_total
            print *, "Total Removal Rate      = ", removed_total
            print *, "Total Incoming Current  = ", incoming_total
            print *, "Total Outgoing Current  = ", escaped_total

        end if
        print *, " "

    end subroutine neutron_balance_diffusion

    subroutine reaction_rate_diffusion(Properties, Results, N)

        type(PropertiesTypeD)   :: Properties
        type(ResultsTypeD)      :: Results
        type(NTypeD)            :: N

        real(kind=8), dimension(N%Material) :: Volume, Mean_Flux, fission_rate, absorption_rate

        integer, dimension(N%Material) :: Number_of_Elements

        integer :: i, j, k

        print *, "Reaction Rates:"
        print *, " "

        do k = 1, N%Group

            Volume = 0.0_8
            fission_rate = 0.0_8
            absorption_rate = 0.0_8
            Mean_Flux = 0.0_8
            Number_of_Elements = 0

            print *, "Group", k
            print *, " "

            print *, "   Material   ", " No. of Elements", "        Volume   ", "                 Mean Flux   ", "              Absorption Rate   ", "        Fission Rate"
            print *, " "

            do j = 1, N%Material
            
                do i = 1, N%Element

                    if (Properties%Elements(i)%Material == j) then

                        fission_rate(j) = fission_rate(j) + (1/Results%k_eff)*Properties%Elements(i)%Sigma_f(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(k,:)))

                        absorption_rate(j) = absorption_rate(j) + Properties%Elements(i)%Sigma_a(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Flux(k,:)))

                        Mean_Flux(j) = Mean_Flux(j) + sum(Properties%Elements(i)%Flux(k,:))*Properties%Elements(i)%Volume/Properties%Elements(i)%Number_of_Nodes

                        Volume(j) = Volume(j) + Properties%Elements(i)%Volume

                        Number_of_Elements(j) = Number_of_Elements(j) + 1

                    end if

                end do

            end do

            if (Properties%g == 1) then

                fission_rate = fission_rate*2.0_8*PI
                absorption_rate = absorption_rate*2.0_8*PI

            else if (Properties%g == 2) then

                fission_rate = fission_rate*4.0_8*PI
                absorption_rate = absorption_rate*4.0_8*PI

            end if

            do j = 1, N%Material

                print *, j, Number_of_Elements(j), '           ',Volume(j), Mean_Flux(j), absorption_rate(j), fission_rate(j)

            end do

            print *, " "

        end do

    end subroutine reaction_rate_diffusion

    subroutine fixed_source_output_diffusion(Properties, N)

        type(PropertiesTypeD)   :: Properties
        type(NTypeD)            :: N

        real(kind=8), dimension(N%Material) :: Volume, Mean_Flux

        integer, dimension(N%Material) :: Number_of_Elements

        integer :: i, j, k

        print *, "Fixed Source Outputs:"
        print *, " "

        do k = 1, N%Group

            Volume = 0.0_8
            Mean_Flux = 0.0_8
            Number_of_Elements = 0

            print *, "Group", k
            print *, " "

            print *, "   Material   ", " No. of Elements", "        Volume   ", "                 Mean Flux   "
            print *, " "

            do j = 1, N%Material
            
                do i = 1, N%Element

                    if (Properties%Elements(i)%Material == j) then

                        Mean_Flux(j) = Mean_Flux(j) + sum(Properties%Elements(i)%Flux(k,:))*Properties%Elements(i)%Volume/Properties%Elements(i)%Number_of_Nodes

                        Volume(j) = Volume(j) + Properties%Elements(i)%Volume

                        Number_of_Elements(j) = Number_of_Elements(j) + 1

                    end if

                end do

            end do

            do j = 1, N%Material

                print *, j, Number_of_Elements(j), '           ',Volume(j), Mean_Flux(j)

            end do

            print *, " "

        end do

    end subroutine fixed_source_output_diffusion

end module