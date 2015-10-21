#!/bin/tcsh
ls results/DNA*/*/*.tab > results/DNA_upto30_tabfiles.list

## Put all of the DNA results tabs together (but skip the header in all but the first).
set started = 0
foreach tabfilename ( "`cat results/DNA_upto30_tabfiles.list`" )
    echo $tabfilename
    if ( $started == 0 ) then
      cat $tabfilename > results/DNA_upto30.tab
      set started = 1
    else
      set inheader = 1
      foreach line ( "`cat $tabfilename`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo $line >> results/DNA_upto30.tab
        endif
      end
    endif
end
# ls results/AA*/*/*.tab > results/AA_upto30_tabfiles.list
    
