#!/bin/tcsh
ls results_revised/AA*/*/*.tab > results_revised/AA_tabfiles.list

foreach tabfilename ( "`cat results_revised/AA_tabfiles.list`" )
    #echo $tabfilename
    perl ProfillicSimulationTabProcessor.pl -sSfzi -O "${tabfilename}.processed.alltogether.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzit -O "${tabfilename}.processed.alltogether.truesseparated.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzi -x 0 -X 0 -O "${tabfilename}.processed.evenstart.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzit -x 0 -X 0 -O "${tabfilename}.processed.evenstart.truesseparated.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzi -x 1 -X 4 -O "${tabfilename}.processed.priorstart.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzit -x 1 -X 4 -O "${tabfilename}.processed.priorstart.truesseparated.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzi -x 5 -X 8 -O "${tabfilename}.processed.uniformstart.out" "$tabfilename"
    perl ProfillicSimulationTabProcessor.pl -sSfzit -x 5 -X 8 -O "${tabfilename}.processed.uniformstart.truesseparated.out" "$tabfilename"
end

## Put all of the AA results tabs together, by tab processing option (but skip the header in all but the first).
set started = 0
foreach tabfilename ( "`cat results_revised/AA_tabfiles.list`" )
    #echo $tabfilename
    if ( $started == 0 ) then
        cat "${tabfilename}.processed.alltogether.out" > "results_revised/AA.processed.alltogether.out"
        cat "${tabfilename}.processed.alltogether.truesseparated.out" > "results_revised/AA.processed.alltogether.truesseparated.out"
        cat "${tabfilename}.processed.evenstart.out" > "results_revised/AA.processed.evenstart.out"
        cat "${tabfilename}.processed.evenstart.truesseparated.out" > "results_revised/AA.processed.evenstart.truesseparated.out"
        cat "${tabfilename}.processed.priorstart.out" > "results_revised/AA.processed.priorstart.out"
        cat "${tabfilename}.processed.priorstart.truesseparated.out" > "results_revised/AA.processed.priorstart.truesseparated.out"
        cat "${tabfilename}.processed.uniformstart.out" > "results_revised/AA.processed.uniformstart.out"
        cat "${tabfilename}.processed.uniformstart.truesseparated.out" > "results_revised/AA.processed.uniformstart.truesseparated.out"
        set started = 1
    else
      ## alltogether
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.alltogether.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.alltogether.out
        endif
      end
      ## alltogether.truesseparated
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.alltogether.truesseparated.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.alltogether.truesseparated.out
        endif
      end

      ## evenstart
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.evenstart.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.evenstart.out
        endif
      end
      ## evenstart.truesseparated
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.evenstart.truesseparated.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.evenstart.truesseparated.out
        endif
      end

      ## priorstart
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.priorstart.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.priorstart.out
        endif
      end
      ## priorstart.truesseparated
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.priorstart.truesseparated.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.priorstart.truesseparated.out
        endif
      end


      ## uniformstart
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.uniformstart.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.uniformstart.out
        endif
      end
      ## uniformstart.truesseparated
      set inheader = 1
      foreach line ( "`cat ${tabfilename}.processed.uniformstart.truesseparated.out`" )
        if ( $inheader > 0 ) then
          set inheader = 0
        else
          echo "$line" >> results_revised/AA.processed.uniformstart.truesseparated.out
        endif
      end
    endif
end

    
