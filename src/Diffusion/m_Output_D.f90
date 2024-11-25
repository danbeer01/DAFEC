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
        real :: origin(N%D)

        top_right_loc = 1


        do j = 1,N%N
            if(N%D == 1) then
                if (this%Nodes(j,1) == Prop%Length(1)) top_right_loc = j
            end if
            if(N%D == 2) then
                if (this%Nodes(j,1) == Prop%Length(1) .and. this%Nodes(j,2) == Prop%Length(2)) top_right_loc = j
            end if
            if(N%D == 3) then
                if (this%Nodes(j,1) == Prop%Length(1) .and. this%Nodes(j,2) == Prop%Length(2) .and. this%Nodes(j,3) == Prop%Length(3)) top_right_loc = j
            end if
        end do

        origin = 0.0

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
        write(*,110) " Position:" , [Prop%Length]
        110 format (A,F8.4,F8.4,F8.4)
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Flux(k,top_right_loc)
        end do
        print *, " "
        print *, "Total ", sum(Results%Flux(:,top_right_loc))
        print *, " "

    end subroutine print_flux_diffusion

    subroutine neutron_balance_diffusion(Properties, Results, N, this)

        type(PropertiesTypeD)   :: Properties
        type(ResultsTypeD)      :: Results
        type(NTypeD)            :: N
        type(Mesh)              :: this

        real(kind=8) :: fission_source = 0.0_8, source_total = 0.0_8
        real(kind=8) :: absorbed = 0.0_8, removed_total = 0.0_8, scattered_out = 0.0_8, scattered_in = 0.0_8
        real(kind=8) :: escaped = 0.0_8, escaped_total = 0.0_8, incoming_total = 0.0_8
        real(kind=8) :: alpha = 0.0_8
        real(kind=8), dimension(N%Group,N%Element) :: Element_Flux
        integer :: i, j, k, Nodes_per_Element

        Element_Flux = 0.0_8

        print *, "Neutron Balance:"
        print *, " "

        if (Properties%LBC == 3 .and. Properties%RBC == 3) alpha = Properties%Alpha

        do k = 1,N%Group

            do i = 1,N%Element
    
                Nodes_per_Element = size(Properties%Elements(i)%Coordinates,1)
    
                do j = 1, Nodes_per_Element
    
                    Element_Flux(k,i) = Element_Flux(k,i) + (Results%Flux(k,this%Cell_Pointers(i,j)))/Nodes_per_Element
    
                end do

            end do
    
        end do

        do k = 1, N%Group
            fission_source = 0.0_8
            absorbed = 0.0_8
            scattered_out = 0.0_8
            scattered_in = 0.0_8

            print *, " Group", k
            print *, " "

            do i = 1, N%Element

                Nodes_per_Element = Properties%Elements(i)%Number_of_Nodes

                do j = 1, N%Group

                    fission_source = fission_source + (1/Results%k_eff)*Properties%Chi(k)*Properties%Elements(i)%Sigma_f(j)*Properties%Elements(i)%Volume*Element_Flux(j,i)

                end do

                do j = k,2,-1 ! Upscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(j-1,k)*Element_Flux(j-1,i)*Properties%Elements(i)%Volume

                end do

                do j = k,N%Group-1 ! Downscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(j+1,k)*Element_Flux(j+1,i)*Properties%Elements(i)%Volume

                end do

                absorbed = absorbed + Properties%Elements(i)%Sigma_a(k)*Properties%Elements(i)%Volume*Element_Flux(k,i)

                do j = 1, N%Group

                    if(j /= k) scattered_out = scattered_out + Properties%Elements(i)%Sigma_s(k,j)*Properties%Elements(i)%Volume*Element_Flux(k,i)

                end do

            end do

            source_total = source_total + fission_source + scattered_in
            print *, "  Source Input      = ", fission_source + scattered_in

            removed_total = removed_total + absorbed + scattered_out
            print *, "  Removal Rate      = ", absorbed + scattered_out

            escaped = fission_source + scattered_in - (absorbed + scattered_out)
            print *, "  Incoming Current  = ", escaped*alpha/(1+alpha)

            if (Properties%LBC == 3) then
                incoming_total = incoming_total + escaped*alpha/(1+alpha)
                escaped = escaped/(1+alpha)
            end if
            
            escaped_total = escaped_total + escaped
            print *, "  Outgoing Current  = ", escaped

            print *, "  Excess Inflow     = ", fission_source + scattered_in - escaped - (absorbed + scattered_out)
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

    subroutine reaction_rate_diffusion(Properties, Results, N, this)

        type(PropertiesTypeD)   :: Properties
        type(ResultsTypeD)      :: Results
        type(NTypeD)            :: N
        type(Mesh)              :: this

        real(kind=8), dimension(N%Material) :: Volume, Mean_Flux, fission_rate, absorption_rate

        real(kind=8) :: Element_Flux

        integer, dimension(N%Material) :: Number_of_Elements

        integer :: i, j, k, a, Nodes_per_Element

        print *, "Reaction Rates:"
        print *, " "

        do k = 1, N%Group

            Volume = 0.0_8
            Element_Flux = 0.0_8
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

                        Nodes_per_Element = Properties%Elements(i)%Number_of_Nodes

                        Element_Flux = 0.0_8

                        do a = 1, Nodes_per_Element
            
                            Element_Flux = Element_Flux + (Results%Flux(k,this%Cell_Pointers(i,a)))/Nodes_per_Element
            
                        end do

                        fission_rate(j) = fission_rate(j) + (1/Results%k_eff)*Properties%Elements(i)%Sigma_f(k)*Properties%Elements(i)%Volume*Element_Flux

                        absorption_rate(j) = absorption_rate(j) + Properties%Elements(i)%Sigma_a(k)*Properties%Elements(i)%Volume*Element_Flux

                        Mean_Flux(j) = Mean_Flux(j) + Element_Flux*Properties%Elements(i)%Volume

                        Volume(j) = Volume(j) + Properties%Elements(i)%Volume

                        Number_of_Elements(j) = Number_of_Elements(j) + 1

                    end if

                end do

                if (Volume(j) > 0.0_8) Mean_Flux(j) = Mean_Flux(j)/Volume(j)

            end do

            do j = 1, N%Material

                print *, j, Number_of_Elements(j), '           ',Volume(j), Mean_Flux(j), absorption_rate(j), fission_rate(j)

            end do

            print *, " "

        end do

    end subroutine reaction_rate_diffusion

end module