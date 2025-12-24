program ramgeo_a100_optimized
c
c     RAMGEO A100 GPU优化版本
c     基于原始ramgeo1.5.f，针对NVIDIA A100 GPU优化
c
c     ******************************************************************
c     ***** Range-dependent Acoustic Model, GPU Optimized Version *****
c     ******************************************************************
c     
      use cudafor
      implicit none
      
      complex ci,ksq,ksqb,u,v,r1,r2,r3,s1,s2,s3,pd1,pd2
      real k0
	complex ksqw
c
c     mr=bathymetry points, mz=depth grid, mp=pade terms.
c
      parameter (mr=5000,mz=999999,mp=10)
      
      ! 使用托管内存（unified memory）
      complex, allocatable, managed :: ksq(:), ksqw(:), ksqb(:)
      complex, allocatable, managed :: u(:), v(:)
      complex, allocatable, managed :: r1(:,:), r2(:,:), r3(:,:)
      complex, allocatable, managed :: s1(:,:), s2(:,:), s3(:,:)
      complex, allocatable, managed :: pd1(:), pd2(:)
      real, allocatable, managed :: tlg(:)
      real, allocatable, managed :: rb(:), zb(:), cw(:), cb(:)
      real, allocatable, managed :: rhob(:), attn(:), alpw(:), alpb(:)
      real, allocatable, managed :: f1(:), f2(:), f3(:)
      
      integer :: mdr, ndr, ndz, iz, nzplt, lz, ib, ir, nz, np, ns
      real :: dir, dr, dz, pi, eta, eps, omega, rmax, c0, r, rp, rs
      integer :: alloc_stat, cuda_stat
      type(cudaDeviceProp) :: prop
      
c
      open(unit=1,status='old',file='ramgeo.in')
      open(unit=2,status='unknown',file='tl.line')
      open(unit=3,status='unknown',file='tl.grid')
c
c     初始化CUDA
      cuda_stat = cudaGetDeviceCount(alloc_stat)
      if (cuda_stat /= cudaSuccess) then
        write(*,*) "警告: 无法获取CUDA设备数量，使用CPU模式"
      else if (alloc_stat > 0) then
        cuda_stat = cudaSetDevice(0)
        cuda_stat = cudaGetDeviceProperties(prop, 0)
        write(*,*) "使用GPU: ", trim(prop%name)
        write(*,*) "GPU内存: ", prop%totalGlobalMem/1024/1024, " MB"
      else
        write(*,*) "警告: 未找到CUDA设备，使用CPU模式"
      endif
c
c     动态分配内存
      allocate(ksq(mz), ksqw(mz), ksqb(mz), stat=alloc_stat)
      allocate(u(mz), v(mz), tlg(mz), stat=alloc_stat)
      allocate(r1(mz,mp), r2(mz,mp), r3(mz,mp), stat=alloc_stat)
      allocate(s1(mz,mp), s2(mz,mp), s3(mz,mp), stat=alloc_stat)
      allocate(pd1(mp), pd2(mp), stat=alloc_stat)
      allocate(rb(mr), zb(mr), cw(mz), cb(mz), stat=alloc_stat)
      allocate(rhob(mz), attn(mz), alpw(mz), alpb(mz), stat=alloc_stat)
      allocate(f1(mz), f2(mz), f3(mz), stat=alloc_stat)
      
      if (alloc_stat /= 0) then
        write(*,*) "错误: 内存分配失败"
        stop
      endif
c
c     使用数据区域管理GPU内存
      !$acc data copyin(ksq, ksqw, ksqb, u, v, r1, r2, r3, s1, s2, s3)
      !$acc& copyin(pd1, pd2, tlg, rb, zb, cw, cb, rhob, attn)
      !$acc& copyin(alpw, alpb, f1, f2, f3)
c
      call setup_gpu(mr, mz, nz, mp, np, ns, mdr, ndr, ndz, iz, 
     >   nzplt, lz, ib, ir, dir, dr, dz, pi, eta, eps, omega, rmax, 
     >   c0, k0, ci, r, rp, rs, rb, zb, cw, cb, rhob, attn, alpw, 
     >   alpb, ksq, ksqw, ksqb, f1, f2, f3, u, v, r1, r2, r3, s1, 
     >   s2, s3, pd1, pd2, tlg)
c
c     March the acoustic field out in range.
c
    1 call updat_gpu(mr, mz, nz, mp, np, iz, ib, dr, dz, eta, 
     >   omega, rmax, c0, k0, ci, r, rp, rs, rb, zb, cw, cb, rhob, 
     >   attn, alpw, alpb, ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, 
     >   s1, s2, s3, pd1, pd2)
     
      call solve_gpu(mz, nz, mp, np, iz, u, v, r1, r2, r3, s1, s2, s3)
      
      r = r + dr
      
      call outpt_gpu(mz, mdr, ndr, ndz, iz, nzplt, lz, ir, dir, 
     >   eps, r, f3, u, tlg, k0)
     
      if (r .lt. rmax) go to 1
c
      !$acc end data
c
      close(1)
      close(2)
      close(3)
c
c     清理内存
      deallocate(ksq, ksqw, ksqb, u, v, tlg)
      deallocate(r1, r2, r3, s1, s2, s3, pd1, pd2)
      deallocate(rb, zb, cw, cb, rhob, attn, alpw, alpb, f1, f2, f3)
c
c     重置CUDA设备
      if (alloc_stat > 0) then
        cuda_stat = cudaDeviceSynchronize()
        cuda_stat = cudaDeviceReset()
      endif
c
      write(*,*) "计算完成"
      stop
      end
c
c     Initialize the parameters, acoustic field, and matrices.
c
      subroutine setup_gpu(mr, mz, nz, mp, np, ns, mdr, ndr, ndz, iz,
     >   nzplt, lz, ib, ir, dir, dr, dz, pi, eta, eps, omega, rmax, 
     >   c0, k0, ci, r, rp, rs, rb, zb, cw, cb, rhob, attn, alpw, 
     >   alpb, ksq, ksqw, ksqb, f1, f2, f3, u, v, r1, r2, r3, s1, s2,
     >   s3, pd1, pd2, tlg)
     
      use cudafor
      complex ci, u(mz), v(mz), ksq(mz), ksqb(mz), r1(mz,mp),
     >   r2(mz,mp), r3(mz,mp), s1(mz,mp), s2(mz,mp), s3(mz,mp), 
     >   pd1(mp), pd2(mp)
      complex ksqw(mz)
      real k0, rb(mr), zb(mr), cw(mz), cb(mz), rhob(mz), attn(mz),
     >   alpw(mz), alpb(mz), f1(mz), f2(mz), f3(mz), tlg(mz)
      integer mdr, ndr, ndz, iz, nzplt, lz, ib, ir, nz, np, ns
      real dir, dr, dz, pi, eta, eps, omega, rmax, c0, r, rp, rs
      real freq, zs, zr, zmax, zmplt, ri, z
      integer i, j, is
c
      read(1,*)
      read(1,*) freq, zs, zr
      read(1,*) rmax, dr, ndr
      read(1,*) zmax, dz, ndz, zmplt
      read(1,*) c0, np, ns, rs
c
      i = 1
    1 read(1,*) rb(i), zb(i)
      if (rb(i) .lt. 0.0) go to 2
      i = i + 1
      go to 1
    2 rb(i) = 2.0 * rmax
      zb(i) = zb(i-1)
c
      pi = 4.0 * atan(1.0)
      ci = cmplx(0.0, 1.0)
      eta = 1.0 / (40.0 * pi * alog10(exp(1.0)))
      eps = 1.0e-20
      ib = 1
      mdr = 0
      r = dr
      omega = 2.0 * pi * freq
      ri = 1.0 + zr / dz
      ir = ifix(ri)
      dir = ri - float(ir)
      k0 = omega / c0
      nz = zmax / dz - 0.5
      nzplt = zmplt / dz - 0.5
      z = zb(1)
      iz = 1.0 + z / dz
      iz = max(2, iz)
      iz = min(nz, iz)
      if (rs .lt. dr) rs = 2.0 * rmax
c
      if (nz + 2 .gt. mz) then
        write(*,*) '需要增加参数 mz 到 ', nz + 2
        stop
      end if
      if (np .gt. mp) then
        write(*,*) '需要增加参数 mp 到 ', np
        stop
      end if
      if (i .gt. mr) then
        write(*,*) '需要增加参数 mr 到 ', i
        stop
      end if
c
      !$acc kernels
      do j = 1, mp
        r3(1,j) = 0.0
        r1(nz+2,j) = 0.0
      end do
      
      do i = 1, nz+2
        u(i) = 0.0
        v(i) = 0.0
      end do
      !$acc end kernels
      
      lz = 0
      do i = ndz, nzplt, ndz
        lz = lz + 1
      end do
      write(3,*) lz
c
c     The initial profiles and starting field.
c
      call profl_gpu(mz, nz, ci, dz, eta, omega, rmax, c0, k0, rp,
     >   cw, cb, rhob, attn, alpw, alpb, ksqw, ksqb)
     
      call selfs_gpu(mz, nz, mp, np, ns, iz, zs, dr, dz, pi, c0, k0,
     >   rhob, alpw, alpb, ksq, ksqw, ksqb, f1, f2, f3, u, v, r1, r2,
     >   r3, s1, s2, s3, pd1, pd2)
     
      call outpt_gpu(mz, mdr, ndr, ndz, iz, nzplt, lz, ir, dir, eps,
     >   r, f3, u, tlg, k0)
c
c     The propagation matrices.
c
      call epade(mp, np, ns, 1, k0, c0, dr, pd1, pd2)
      call matrc_gpu(mz, nz, mp, np, iz, iz, dz, k0, rhob, alpw, alpb,
     >   ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
c
      return
      end
c
c     Set up the profiles.
c
      subroutine profl_gpu(mz, nz, ci, dz, eta, omega, rmax, c0, k0, rp,
     >   cw, cb, rhob, attn, alpw, alpb, ksqw, ksqb)
     
      use cudafor
      complex ci, ksqb(mz), ksqw(mz)
      real k0, cw(mz), cb(mz), rhob(mz), attn(mz), alpw(mz), alpb(mz)
      real, allocatable, managed :: cwa(:)
      integer i
c
      call zread_gpu(mz, nz, dz, cw)
      allocate(cwa(mz), stat=i)
      if (i /= 0) then
        write(*,*) "错误: 无法分配cwa数组"
        stop
      endif
      
      call zread_gpu(mz, nz, dz, cwa)
      call zread_gpu(mz, nz, dz, cb)
      call zread_gpu(mz, nz, dz, rhob)
      call zread_gpu(mz, nz, dz, attn)
      
      rp = 2.0 * rmax
      read(1,*,end=1) rp
c
    1 continue
      
      !$acc kernels present(cw, cb, rhob, attn, alpw, alpb, ksqw, ksqb)
      !$acc& present(cwa)
      do i = 1, nz+2
        ksqw(i) = ((omega / cw(i)) * (1.0 + ci * eta * cwa(i)))**2 - k0**2
        ksqb(i) = ((omega / cb(i)) * (1.0 + ci * eta * attn(i)))**2 - k0**2
        alpw(i) = sqrt(cw(i) / c0)
        alpb(i) = sqrt(rhob(i) * cb(i) / c0)
      end do
      !$acc end kernels
c
      deallocate(cwa)
      
      return
      end
c
c     Profile reader and interpolator.
c
      subroutine zread_gpu(mz, nz, dz, prof)
      real prof(mz)
      real zi, profi
      integer i, j, k, iold
c
      !$acc kernels
      do i = 1, nz+2
        prof(i) = -1.0
      end do
      !$acc end kernels
      
      read(1,*) zi, profi
      prof(1) = profi
      i = 1.5 + zi / dz
      prof(i) = profi
      iold = i
    2 read(1,*) zi, profi
      if (zi .lt. 0.0) go to 3
      i = 1.5 + zi / dz
      if (i .eq. iold) i = i + 1
      prof(i) = profi
      iold = i
      go to 2
    3 prof(nz+2) = prof(i)
      i = 1
      j = 1
    4 i = i + 1
      if (prof(i) .lt. 0.0) go to 4
      if (i - j .eq. 1) go to 6
      
      !$acc kernels present(prof)
      do k = j+1, i-1
        prof(k) = prof(j) + float(k-j) * (prof(i) - prof(j)) / float(i-j)
      end do
      !$acc end kernels
      
    6 j = i
      if (j .lt. nz+2) go to 4
c
      return
      end
c
c     The tridiagonal matrices.
c
      subroutine matrc_gpu(mz, nz, mp, np, iz, jz, dz, k0, rhob, alpw,
     >   alpb, ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3,
     >   pd1, pd2)
     
      use cudafor
      complex d1, d2, d3, rfact, ksq(mz), ksqb(mz), r1(mz,mp),
     >   r2(mz,mp), r3(mz,mp), s1(mz,mp), s2(mz,mp), s3(mz,mp),
     >   pd1(mp), pd2(mp)
      real k0, rhob(mz), f1(mz), f2(mz), f3(mz), alpw(mz), alpb(mz)
      complex ksqw(mz)
      real a1, a2, a3, cfact, dfact, c1, c2, c3
      integer i, j, ii
c
      a1 = k0**2 / 6.0
      a2 = 2.0 * k0**2 / 3.0
      a3 = k0**2 / 6.0
      cfact = 0.5 / dz**2
      dfact = 1.0 / 12.0
c
      !$acc kernels present(f1, f2, f3, ksq, alpw, alpb, rhob, ksqw, ksqb)
      do i = 1, iz
        f1(i) = 1.0 / alpw(i)
        f2(i) = 1.0
        f3(i) = alpw(i)
        ksq(i) = ksqw(i)
      end do
      
      ii = 1
      do i = iz+1, nz+2
        f1(i) = rhob(ii) / alpb(ii)
        f2(i) = 1.0 / rhob(ii)
        f3(i) = alpb(ii)
        ksq(i) = ksqb(ii)
        ii = ii + 1
      end do
      !$acc end kernels
c
      !$acc parallel loop collapse(2) present(f1, f2, f3, ksq, r1, r2, r3)
      !$acc& present(s1, s2, s3, pd1, pd2)
      !$acc& private(c1, c2, c3, d1, d2, d3)
      do i = 2, nz+1
        do j = 1, np
c
c         Discretization by Galerkin's method.
c
          c1 = cfact * f1(i) * (f2(i-1) + f2(i)) * f3(i-1)
          c2 = -cfact * f1(i) * (f2(i-1) + 2.0*f2(i) + f2(i+1)) * f3(i)
          c3 = cfact * f1(i) * (f2(i) + f2(i+1)) * f3(i+1)
          
          d1 = c1 + dfact * (ksq(i-1) + ksq(i))
          d2 = c2 + dfact * (ksq(i-1) + 6.0*ksq(i) + ksq(i+1))
          d3 = c3 + dfact * (ksq(i) + ksq(i+1))
c
          r1(i,j) = a1 + pd2(j) * d1
          r2(i,j) = a2 + pd2(j) * d2
          r3(i,j) = a3 + pd2(j) * d3
          s1(i,j) = a1 + pd1(j) * d1
          s2(i,j) = a2 + pd1(j) * d2
          s3(i,j) = a3 + pd1(j) * d3
        end do
      end do
      !$acc end parallel
c
c     The matrix decomposition.
c
      !$acc parallel loop collapse(2) present(r1, r2, r3, s1, s2, s3)
      !$acc& private(rfact)
      do j = 1, np
        do i = 2, nz+1
          rfact = 1.0 / (r2(i,j) - r1(i,j) * r3(i-1,j))
          r1(i,j) = r1(i,j) * rfact
          r3(i,j) = r3(i,j) * rfact
          s1(i,j) = s1(i,j) * rfact
          s2(i,j) = s2(i,j) * rfact
          s3(i,j) = s3(i,j) * rfact
        end do
      end do
      !$acc end parallel
c
      return
      end
c
c     The tridiagonal solver.
c
      subroutine solve_gpu(mz, nz, mp, np, iz, u, v, r1, r2, r3, s1, s2, s3)
      use cudafor
      complex u(mz), v(mz), r1(mz,mp), r2(mz,mp), r3(mz,mp),
     >   s1(mz,mp), s2(mz,mp), s3(mz,mp)
      complex, parameter :: eps_cmplx = (1.0e-30, 0.0)
      integer j, i
c
      do j = 1, np
c
        !$acc kernels present(u, v, s1, s2, s3)
        do i = 2, nz+1
          v(i) = s1(i,j) * u(i-1) + s2(i,j) * u(i) + s3(i,j) * u(i+1) + eps_cmplx
        end do
        !$acc end kernels
c
        !$acc kernels present(v, r1)
        do i = 3, nz+1
          v(i) = v(i) - r1(i,j) * v(i-1) + eps_cmplx
        end do
        !$acc end kernels
c
        !$acc kernels present(u, v, r3)
        do i = nz, 2, -1
          u(i) = v(i) - r3(i,j) * u(i+1) + eps_cmplx
        end do
        !$acc end kernels
        
      end do
c
      return
      end
c
c     Matrix updates.
c
      subroutine updat_gpu(mr, mz, nz, mp, np, iz, ib, dr, dz, eta, omega,
     >   rmax, c0, k0, ci, r, rp, rs, rb, zb, cw, cb, rhob, attn, alpw,
     >   alpb, ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
     
      use cudafor
      complex ci, ksq(mz), ksqb(mz), r1(mz,mp), r2(mz,mp), r3(mz,mp),
     >   s1(mz,mp), s2(mz,mp), s3(mz,mp), pd1(mp), pd2(mp)
      real k0, rb(mr), zb(mr), attn(mz), cb(mz), rhob(mz), cw(mz),
     >   f1(mz), f2(mz), f3(mz), alpb(mz), alpw(mz)
      complex ksqw(mz)
      real z
      integer jz
c
c     Varying bathymetry.
c
      if (r .ge. rb(ib+1)) ib = ib + 1
      jz = iz
      z = zb(ib) + (r + 0.5*dr - rb(ib)) * (zb(ib+1) - zb(ib)) /
     >   (rb(ib+1) - rb(ib))
      iz = 1.0 + z / dz
      iz = max(2, iz)
      iz = min(nz, iz)
      
      if (iz .ne. jz) then
        call matrc_gpu(mz, nz, mp, np, iz, jz, dz, k0, rhob, alpw, alpb,
     >     ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
      endif
c
c     Varying profiles.
c
      if (r .ge. rp) then
        call profl_gpu(mz, nz, ci, dz, eta, omega, rmax, c0, k0, rp,
     >     cw, cb, rhob, attn, alpw, alpb, ksqw, ksqb)
        call matrc_gpu(mz, nz, mp, np, iz, iz, dz, k0, rhob, alpw, alpb,
     >     ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
      endif
c
c     Turn off the stability constraints.
c
      if (r .ge. rs) then
        ns = 0
        rs = 2.0 * rmax
        call epade(mp, np, ns, 1, k0, c0, dr, pd1, pd2)
        call matrc_gpu(mz, nz, mp, np, iz, iz, dz, k0, rhob, alpw, alpb,
     >     ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
      endif
c
      return
      end
c
c     The self-starter.
c
      subroutine selfs_gpu(mz, nz, mp, np, ns, iz, zs, dr, dz, pi, c0, k0,
     >   rhob, alpw, alpb, ksq, ksqw, ksqb, f1, f2, f3, u, v, r1, r2, r3,
     >   s1, s2, s3, pd1, pd2)
     
      use cudafor
      complex u(mz), v(mz), ksq(mz), ksqb(mz), r1(mz,mp), r2(mz,mp),
     >   r3(mz,mp), s1(mz,mp), s2(mz,mp), s3(mz,mp), pd1(mp), pd2(mp)
      real k0, rhob(mz), alpw(mz), alpb(mz), f1(mz), f2(mz), f3(mz)
      complex ksqw(mz)
      real si, dis
      integer is
c
c     Conditions for the delta function.
c
      si = 1.0 + zs / dz
      is = ifix(si)
      dis = si - float(is)
      
      !$acc kernels present(u, alpw)
      u(is) = (1.0 - dis) * sqrt(2.0 * pi / k0) / (dz * alpw(is))
      u(is+1) = dis * sqrt(2.0 * pi / k0) / (dz * alpw(is))
      !$acc end kernels
c
c     Divide the delta function by (1-X)**2 to get a smooth rhs.
c
      pd1(1) = 0.0
      pd2(1) = -1.0
      call matrc_gpu(mz, nz, mp, 1, iz, iz, dz, k0, rhob, alpw, alpb,
     >   ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
      
      call solve_gpu(mz, nz, mp, 1, iz, u, v, r1, r2, r3, s1, s2, s3)
      call solve_gpu(mz, nz, mp, 1, iz, u, v, r1, r2, r3, s1, s2, s3)
c
c     Apply the operator (1-X)**2*(1+X)**(-1/4)*exp(ci*k0*r*sqrt(1+X)).
c
      call epade(mp, np, ns, 2, k0, c0, dr, pd1, pd2)
      call matrc_gpu(mz, nz, mp, np, iz, iz, dz, k0, rhob, alpw, alpb,
     >   ksq, ksqw, ksqb, f1, f2, f3, r1, r2, r3, s1, s2, s3, pd1, pd2)
      
      call solve_gpu(mz, nz, mp, np, iz, u, v, r1, r2, r3, s1, s2, s3)
c
      return
      end
c
c     Output transmission loss.
c
      subroutine outpt_gpu(mz, mdr, ndr, ndz, iz, nzplt, lz, ir, dir,
     >   eps, r, f3, u, tlg, k0)
     
      use cudafor
      complex ur, u(mz), ci
      complex, allocatable, managed :: tlgg(:)
      real f3(mz), tlg(mz), k0
      real pi, tl
      integer j, i, alloc_stat
c
      pi = 4.0 * atan(1.0)
      ci = cmplx(0.0, 1.0)
      
      ! 使用托管内存分配
      allocate(tlgg(mz), stat=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "警告: 无法分配tlgg内存，跳过此输出"
        return
      endif
c
      !$acc kernels present(u, f3)
      ur = (1.0 - dir) * f3(ir) * u(ir) + dir * f3(ir+1) * u(ir+1)
      !$acc end kernels
      
      tl = -20.0 * alog10(cabs(ur) + eps) + 10.0 * alog10(r + eps)
      
      write(2,*) r, tl, -real(exp(-pi*ci/4+ci*k0*r)*ci*ur/sqrt(r+eps)),
     &   imag(exp(-pi*ci/4+ci*k0*r)*ci*ur/sqrt(r+eps))
c
      mdr = mdr + 1
      if (mdr .eq. ndr) then
        mdr = 0
        j = 0
        
        !$acc parallel loop present(u, f3, tlg, tlgg) private(ur)
        do i = ndz, nzplt, ndz
          ur = u(i) * f3(i)
          j = j + 1
          tlg(j) = -20.0 * alog10(cabs(ur) + eps) + 10.0 * alog10(r + eps)
          tlgg(j) = -conjg(exp(-pi*ci/4 + ci*k0*r) * ci * ur) / sqrt(r + eps)
        end do
        !$acc end parallel
        
        write(3,*) (tlg(j), real(tlgg(j)), imag(tlgg(j)), j = 1, lz)
      endif
c
      if (allocated(tlgg)) then
        deallocate(tlgg)
      endif
      
      return
      end