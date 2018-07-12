program test_parse
   use parser
   implicit none
   integer:: err, keyval
   character(len=60):: key
   character(len=256):: line

   call getarg(1,line)
   open(unit=10,file=trim(line),status='OLD')
   err = 0
   ! first we have a section header called 'ABBC'
   call parse(10, key, line, err)

   if (trim(key).eq.'SECTION') then
      if (trim(line).eq.'ABBC') then
         Write(*,*) 'Pass.'
      else
         Write(*,*) 'Wrong section name: ', trim(line)
         Write(*,*) 'Fail.'
         stop
      endif
   else
      Write(*,*) 'Wrong entry type: ', trim(key)
      Write(*,*) 'Fail.'
      stop
   endif

   ! Next we have a comment and an empty line which will be discarded.
   ! and after that a key, called 'key1' which should have the value 1234

   call parse(10, key, line, err)
   if (trim(key).eq.'key1') then
      call read_val(line, keyval)
      if (keyval.eq.1234) then
         Write(*,*) 'Pass.'
      else
         Write(*,*) 'Wrong key value: ', keyval
         Write(*,*) 'Fail.'
         stop
      endif
   else
      Write(*,*) 'Wrong key name: ', trim(key)
      Write(*,*) 'Fail.'
      stop
   endif
      

end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
