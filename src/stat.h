c     %--------------------------------%
c     | See stat.doc for documentation |
c     %--------------------------------%
c
c\SCCS Information: @(#) 
c FILE: stat.h   SID: 2.2   DATE OF SID: 11/16/95   RELEASE: 2 
c
c      real       t0, t1, t2, t3, t4, t5
c      save       t0, t1, t2, t3, t4, t5
c
      integer    nopx, nbx, nrorth
cp         , nitref, nrstrt
ct      real       tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
ct     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
ct     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
ct     &           tmvopx, tmvbx, tgetv0, titref, trvec
      common /timing/ 
     &           nopx, nbx, nrorth
cp  , nitref, nrstrt
ct     &           tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv,
ct     &           tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv,
ct     &           tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv,
ct     &           tmvopx, tmvbx, tgetv0, titref, trvec
