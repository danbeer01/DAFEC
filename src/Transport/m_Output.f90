module m_output
!
! Purpose:
! This module contains the subroutines for printing the flux, k_eff and neutron balance
!
! Record of revisions:
! Date          Programmer      Description of change
! ==========    ===========     =====================
! 31/10/2024    Daniel Beer     Original code
!
use m_Read_Properties
use m_Constants_mod
use m_Results
use m_Calculate_mu_w
use m_VTK_Reader

implicit none

contains
    
    subroutine print_flux(Prop, Results, N, this)

        type(PropertiesType)   :: Prop
        type(ResultsType)      :: Results
        type(NType)            :: N
        type(Mesh)             :: this

        integer :: j, k, top_right_loc
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
        print *, k, Results%Scalar_Flux(k,1)
        end do
        print *, " "
        print *, "Total ", sum(Results%Scalar_Flux(:,1))
        print *, " "
        write(*,110) " Position:", this%Nodes((MAXLOC(Results%Scalar_Flux(1,:))),:)
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Scalar_Flux(k,(MAXLOC(Results%Scalar_Flux(1,:))))
        end do
        print *, " "
        print *, "Total ", sum(Results%Scalar_Flux(:,(MAXLOC(Results%Scalar_Flux(1,:)))))
        print *, " "
        write(*,110) " Position:" , [top_right_coord]
        110 format (A,F8.4,F8.4,F8.4)
        print *, " "
        do k = 1, N%Group
        print *, k, Results%Scalar_Flux(k,top_right_loc)
        end do
        print *, " "
        print *, "Total ", sum(Results%Scalar_Flux(:,top_right_loc))
        print *, " "

    end subroutine print_flux

    subroutine neutron_balance(Properties, Results, N)

        type(PropertiesType)   :: Properties
        type(ResultsType)      :: Results
        type(NType)            :: N

        real(kind=8) :: fission_source = 0.0_8, source_total = 0.0_8
        real(kind=8) :: absorbed = 0.0_8, removed_total = 0.0_8, scattered_out = 0.0_8, scattered_in = 0.0_8
        real(kind=8) :: escaped, escaped_total = 0.0_8, incoming_total = 0.0_8
        real(kind=8) :: alpha = 0.0_8
        real(kind=8), dimension(N%Ordinates) :: escaped_ang
        integer :: i, j, k, ang, side_index, neighbour_element, w_index
        real(kind=8), dimension(:), allocatable :: w

        if (N%D == 1) allocate(w(N%Ordinates))
        if (N%D == 2) allocate(w(N%Ordinates/4))
        if (N%D == 3) allocate(w(N%Ordinates/8))

        if (N%D == 1) call Calculate_w_1D(w)
        if (N%D == 2) call Calculate_w(w)
        if (N%D == 3) call Calculate_w(w)

        print *, "Neutron Balance:"
        print *, " "

        if (Properties%LBC == 3) alpha = Properties%Alpha

        do k = 1, N%Group

            fission_source = 0.0_8
            absorbed = 0.0_8
            scattered_out = 0.0_8
            scattered_in = 0.0_8
            escaped_ang = 0.0_8
            escaped = 0.0_8

            print *, " Group", k
            print *, " "

            do i = 1, N%Element

                do j = 1, N%Group

                    fission_source = fission_source + (1/Results%k_eff)*Properties%Chi(k)*Properties%Elements(i)%Sigma_f(j)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(j,:)))

                end do

                do j = k,2,-1 ! Upscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(0,j-1,k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(j-1,:)))

                end do

                do j = k,N%Group-1 ! Downscatter

                    scattered_in = scattered_in + Properties%Elements(i)%Sigma_s(0,j+1,k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(j+1,:)))

                end do

                absorbed = absorbed + Properties%Elements(i)%Sigma_a(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(k,:)))

                do j = 1, N%Group

                    if(j /= k) scattered_out = scattered_out + Properties%Elements(i)%Sigma_s(0,k,j)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(k,:)))

                end do

                do ang = 1, N%Ordinates

                    do side_index = 1, Properties%Elements(i)%Number_of_Sides

                        neighbour_element = Properties%Elements(i)%Neighbours(side_index,1)

                        if (neighbour_element == 0) then

                            escaped_ang(ang) = escaped_ang(ang) + sum(matmul((-Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:)),Properties%Elements(i)%Sides(side_index)%Boundary(k,ang,:)))

                        else

                            escaped_ang(ang) = escaped_ang(ang) + sum(matmul((-Properties%Elements(i)%Sides(side_index)%F_in_Matrix(ang,:,:)),Properties%Elements(neighbour_element)%Flux(k,ang,:)))

                        end if

                    end do

                    escaped_ang(ang) = escaped_ang(ang) + sum(matmul((Properties%Elements(i)%F_out_Matrix(ang,:,:)),Properties%Elements(i)%Flux(k,ang,:)))

                end do

            end do

            do ang = 1, N%Ordinates

                if (N%D == 1) w_index = ang
                if (N%D == 2) w_index = MOD(ang,N%Ordinates/4)
                if (N%D == 3) w_index = MOD(ang,N%Ordinates/8)

                if (w_index == 0) then

                    if (N%D == 2) w_index = N%Ordinates/4
                    if (N%D == 3) w_index = N%Ordinates/8

                end if

                escaped = escaped + (1.0_8/(2.0_8**N%D))*w(w_index)*escaped_ang(ang)

            end do

            if (Properties%g == 1) then

                fission_source = fission_source*2.0_8*PI
                scattered_in = scattered_in*2.0_8*PI
                scattered_out = scattered_out*2.0_8*PI
                absorbed = absorbed*2.0_8*PI
                escaped = escaped*2.0_8*PI

            else if (Properties%g == 2) then

                fission_source = fission_source*4.0_8*PI
                scattered_in = scattered_in*4.0_8*PI
                scattered_out = scattered_out*4.0_8*PI
                absorbed = absorbed*4.0_8*PI
                escaped = escaped*4.0_8*PI

            end if

            source_total = source_total + fission_source + scattered_in
            print *, "  Source Input      = ", fission_source + scattered_in

            removed_total = removed_total + absorbed + scattered_out
            print *, "  Removal Rate      = ", absorbed + scattered_out

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

    end subroutine neutron_balance

    subroutine reaction_rate(Properties, Results, N)

        type(PropertiesType)   :: Properties
        type(ResultsType)      :: Results
        type(NType)            :: N

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

            print *, "   Material   ", "   No. of Elements", "      Volume   ", "                 Mean Flux   ", "              Absorption Rate   ", "        Fission Rate"
            print *, " "

            do j = 1, N%Material
            
                do i = 1, N%Element
                    
                    if (Properties%Elements(i)%Material == j) then

                        fission_rate(j) = fission_rate(j) + (1/Results%k_eff)*Properties%Elements(i)%Sigma_f(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(k,:)))

                        absorption_rate(j) = absorption_rate(j) + Properties%Elements(i)%Sigma_a(k)*sum(matmul(Properties%Elements(i)%A_Matrix,Properties%Elements(i)%Scalar_Flux(k,:)))

                        Mean_Flux(j) = Mean_Flux(j) + sum(Properties%Elements(i)%Scalar_Flux(k,:))*Properties%Elements(i)%Volume/Properties%Elements(i)%Number_of_Nodes

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

    end subroutine reaction_rate

end module m_output