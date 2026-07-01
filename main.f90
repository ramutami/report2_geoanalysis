
module parameters_and_variables_for_simulation
    implicit none
    real(8), parameter :: pi = atan(1.0d0)*4.0d0
    !gloval_parameters
        !global定数：すべてのmoduleのsubroutineで参照可能な変数

        !------------change allowed--------------!
        real(8),parameter :: g = 0.02d0 ![m/s^2]
        real(8),parameter :: H_par = 750.0d0 ![m]
        real(8),parameter :: A = 1.0d3 ![m^2/s]
        real(8),parameter :: Omega = 7.2921159d-5 ![rad/sec]

        real(8),parameter :: delta_t_day = 0.005d0 ![day]
        integer,parameter :: outstep_interval = 100 
        integer,parameter :: timestep_max = 200*365
        integer,parameter :: smoothify_step = 3 !平滑化間隔

        real(8),parameter :: gamma = 0.005d0

        integer,parameter :: width_deg = 60 ![計算空間の横幅(経度,deg)]
        integer,parameter :: minimum_latitude = 20 ![緯度の下限]
        integer,parameter :: height_deg = 20 ![計算空間の縦幅(緯度,deg)]
        integer,parameter :: grid_perdeg = 8 ![1degをgrid_degこの格子に切る]

        real(8),parameter :: sigma_gaussian_x = 2.0d0 !初期条件のガウス分布のx方向のsigma[deg]。
        real(8),parameter :: sigma_gaussian_y = 2.0d0 !初期条件のガウス分布のx方向のsigma[deg]。
        real(8),parameter :: scale_gaussian = 50.0d-2 !初期分布の高さスケール（うまく調整する。）


        !------------change allowed--------------!

        real(8),parameter :: Rearth = 6378137.0d0 ![地球半径 m]
        real(8),parameter :: dy = Rearth*pi/180.0d0/real(grid_perdeg,8) ![y方向のグリッド幅(m)]
        real(8),parameter :: dt = delta_t_day * 86400.0d0 !dt[sec]

        integer,parameter :: imax = width_deg*grid_perdeg+1 
        integer,parameter :: jmax = height_deg*grid_perdeg+1
        real(8),parameter :: width = pi*width_deg/180.0d0 ![計算空間の横幅、rad]
        real(8),parameter :: height = pi*height_deg/180.0d0 ![計算空間の縦幅、rad]

    !

    !global_variables
        !global変数：すべてのmoduleのsubroutineで参照可能な変数

        real(8),allocatable :: un(:,:),up(:,:),u(:,:)
        real(8),allocatable :: vn(:,:),vp(:,:),v(:,:)
        real(8),allocatable :: hn(:,:),hp(:,:),h(:,:)
        real(8),allocatable :: phi(:),lambda(:) !緯度経度[rad]
        real(8),allocatable :: dx(:) !緯度ごとのdx

    !

    contains
    subroutine allocate_variables() 
        implicit none   
        integer :: i,j
        real(8) :: latitude_j,longtitude_i

        allocate(dx(1:jmax))
        allocate(phi(1:jmax))
        allocate(lambda(1:imax))
        allocate(un(1:imax+1,1:jmax),up(1:imax+1,1:jmax),u(1:imax+1,1:jmax))
        allocate(vn(1:imax,1:jmax+1),vp(1:imax,1:jmax+1),v(1:imax,1:jmax+1))
        allocate(hn(1:imax,1:jmax),hp(1:imax,1:jmax),h(1:imax,1:jmax))

        un = 0.0d0
        up = 0.0d0
        u = 0.0d0

        vn = 0.0d0
        vp = 0.0d0
        v = 0.0d0

        hp = 0.0d0
        hn = 0.0d0
        h = 0.0d0
        
        do j = 1,jmax
            latitude_j = (real(minimum_latitude,8)+real(j-1,8)/real(grid_perdeg,8))*pi/180.0d0 !jでの緯度[rad]
            phi(j) = latitude_j
            dx(j) = pi/real(grid_perdeg,8)/180.0d0 * Rearth * cos(latitude_j)
        end do
        
        do i = 1,imax
            longtitude_i = (0.0d0 + real(i-1,8)/real(grid_perdeg,8))*pi/180.0d0
            lambda(i) = longtitude_i
        end do

    end subroutine allocate_variables

end module parameters_and_variables_for_simulation

module initialization
    use parameters_and_variables_for_simulation
    implicit none
    contains

    subroutine initialize_gaussian
        implicit none
        integer :: i_center_gaussian,j_center_gaussian
        real(8) :: haversine_dist 
        real(8) :: val1,val2
        real(8) :: sigma_x2,sigma_y2
        integer :: i,j

        i_center_gaussian = imax
        j_center_gaussian = floor(real(1+jmax,8)/2.0d0)

        sigma_x2 = (sigma_gaussian_x*pi/180.0d0)**2 !sigma_x[rad]^2
        sigma_y2 = (sigma_gaussian_y*pi/180.0d0)**2 !sigma_y[rad]^2
        

        do j = 1,jmax
        do i = 1,imax
            val1 = (lambda(i)-lambda(i_center_gaussian))**2/sigma_x2
            val2 = (phi(j)-phi(j_center_gaussian))**2/sigma_y2
            h(i,j) = scale_gaussian*exp(-0.5d0*(val1+val2))/sqrt(4.0d0*pi*pi*sigma_x2 * sigma_y2)
        end do
        end do

        hp = h
        
    end subroutine initialize_gaussian
end module initialization

module time_integration
    use parameters_and_variables_for_simulation
    implicit none
    contains

    subroutine calculate_u
        implicit none 
        integer :: i,j
        real(8) :: term1,term2,term3

        do j = 2,jmax-1
        do i = 2,imax
            term1 = 2.0d0*Omega*sin(phi(j))*0.25d0*(v(i-1,j)+v(i-1,j+1)+v(i,j)+v(i,j+1))
            term2 = -g*(h(i,j)-h(i-1,j))/dx(j)
            term3 = A*((up(i+1,j)-2.0d0*up(i,j)+up(i-1,j))/dx(j)**2 + (up(i,j+1)-2.0d0*up(i,j)+up(i,j-1))/dy**2)
            un(i,j) = up(i,j)+2.0d0*dt*(term1+term2+term3)
        end do
        end do

        do j=1,jmax
            un(1,j) = 0.0d0
            un(imax+1,j) = 0.0d0
        end do

        do i=1,imax+1
            un(i,1) = 0.0d0
            un(i,jmax) = 0.0d0
        end do

    end subroutine calculate_u

    subroutine calculate_v
        implicit none
        integer :: i,j
        real(8) :: term1,term2,term3

        do j = 2,jmax
        do i = 2,imax-1
            term1 = -2.0d0*Omega*sin(0.5d0*(phi(j-1)+phi(j)))*0.25d0*(u(i,j)+u(i+1,j)+u(i+1,j-1)+u(i,j-1))
            term2 = -g*(h(i,j)-h(i,j-1))/dy
            term3 = A*((vp(i+1,j)-2.0d0*vp(i,j)+vp(i-1,j))/(0.5d0*(dx(j)+dx(j-1)))**2 + (vp(i,j+1)-2.0d0*vp(i,j)+vp(i,j-1))/dy**2 )
            vn(i,j) = vp(i,j)+2.0d0*dt*(term1+term2+term3)   
        end do
        end do

        do j=1,jmax+1
            vn(1,j)=0.0d0
            vn(imax,j)=0.0d0
        end do

        do i=1,imax
            vn(i,1)=0.0d0
            vn(i,jmax+1)=0.0d0
        end do

    end subroutine calculate_v

    subroutine calculate_h
        implicit none
        integer :: i,j

        do j=1,jmax
        do i=1,imax
            hn(i,j) = hp(i,j)+2.0d0*dt*(-H_par*((u(i+1,j)-u(i,j))/dx(j) + (v(i,j+1)-v(i,j))/dy))
        end do
        end do

    end subroutine calculate_h

    subroutine smoothify
        implicit none
        integer :: i,j

        u = u + 0.5d0*gamma*(un-2.0d0*u+up)
        v = v + 0.5d0*gamma*(vn-2.0d0*v+vp)
        h = h + 0.5d0*gamma*(hn-2.0d0*h+hp)

        do j=1,jmax
            u(1,j) = 0.0d0
            u(imax+1,j) = 0.0d0
        end do

        do i=1,imax+1
            u(i,1) = 0.0d0
            u(i,jmax) = 0.0d0
        end do

        do j=1,jmax+1
            v(1,j)=0.0d0
            v(imax,j)=0.0d0
        end do

        do i=1,imax
            v(i,1)=0.0d0
            v(i,jmax+1)=0.0d0
        end do

    end subroutine smoothify

    subroutine move_one_step
        implicit none

        up = u
        u = un
        vp = v
        v = vn
        hp=h
        h=hn

    end subroutine move_one_step

end module time_integration

module output_module
    use parameters_and_variables_for_simulation
    implicit none
    contains

    subroutine writeout_data(time,output_step,timestep)
        implicit none
        real(8),intent(in) :: time
        integer,intent(in) :: output_step,timestep
        character(len=128) :: filename
        integer :: i,j
        real(8) :: longtitude_deg,latitude_deg
        write(filename,fmt='("./dataout/dataout",i6.6,".dat")') output_step
        open(10,file=filename,status='replace',action='write')
        write(10,*) "# time[dat]",time
        write(10,*) "# calculation step",timestep
        write(10,*) "# longtitude[deg]/latitude[deg],h[m]/u[m/s]/v[m/s]"
        !セル中央での値を出力する！
        do i = 1,imax
        do j = 1,jmax
            longtitude_deg = lambda(i)*180.0d0/pi
            latitude_deg = phi(j)*180.0d0/pi
            write(10,*) longtitude_deg,latitude_deg,h(i,j),0.5d0*(u(i,j)+u(i+1,j)),0.5d0*(v(i,j)+v(i,j+1))
        end do
        write(10,*)
        end do
        close(10)
    end subroutine writeout_data

end module output_module

program main
    use parameters_and_variables_for_simulation
    use initialization
    use time_integration
    use output_module
    implicit none 
    integer :: timestep,outstep
    real(8) :: time

    !初期状態
    call allocate_variables()
    call initialize_gaussian()

    outstep = 0
    do timestep = 0,timestep_max-1

        !出力
        if (mod(timestep,outstep_interval)==0) then
            outstep = outstep + 1
            write(*,'(f4.1,A,i4)') real(timestep,8)*100.0/real(timestep_max-1,8),"%, outputstep=",outstep
            time = delta_t_day  * (timestep)
            call writeout_data(time,outstep,timestep)
        end if

        !計算
        call calculate_u()
        call calculate_v()
        call calculate_h()

        !平滑化
        if (mod(timestep,smoothify_step)==0 .and. timestep /= 0) then
            call smoothify()
        end if

        !更新
        call move_one_step()

    end do

end program main
