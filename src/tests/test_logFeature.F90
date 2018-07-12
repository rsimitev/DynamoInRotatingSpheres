! Copyright  L. Silva (lacsilva@gmail.com), 2015
program testLogFeature
   use drs_params

   call logFeature('New feature', 4, 'A new feature I just made up.')
   call logFeature('Another feature', 5, 'A second new feature I just made up.')
   call logFeature('Very looonnnnnnngggg feature', 20, 'A third new feature with a long short name.')
end program
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
