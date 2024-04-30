subroutine readparam(iconf,conf)

   integer iconf(7)
   real*8  conf(2)
   integer ios

   open(unit=5,file='param',iostat=ios,status='old')
   if(ios.ne.0) then
    print *,'Reading error param file - IOSTAT = ',ios
    stop
   endif

   read(5,*)
   read(5,*) iconf(1)
   read(5,*)
   read(5,*) iconf(2)
   read(5,*)
   read(5,*) iconf(3)
   read(5,*)
   read(5,*) iconf(4)
   read(5,*)
   read(5,*) iconf(5)
   read(5,*)
   read(5,*) iconf(6)
   read(5,*)
   read(5,*) iconf(7)
   read(5,*)
   read(5,*) conf(1)
   read(5,*)
   read(5,*) conf(2)

   close(5)

   return
end
