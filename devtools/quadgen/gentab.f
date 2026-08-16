c
c       Dump the chunkmatc auxiliary Galerkin quadrature tables
c       (legeexps_log_lr aux targets + hqsuppauxquad per-target
c       log-singular source rules) as plain text, for formatting into
c       chunkie's chnk.quadgalerkin tables by mktab.py.
c
c       See README.md for build and usage.
c
c       Orders to dump are read from stdin, one integer per line.
c
c       Output format (all nodes/weights mapped to [-1,1]):
c
c         ORDER <n> <naux>
c         AUX
c         <ts_aux(j)> <ws_aux(j)>        (naux lines)
c         SRC <inode> <npts>
c         <xs0(i)> <ws0(i)>              (npts lines)
c         ... (naux SRC blocks)
c
        program gentab
        implicit real *8 (a-h,o-z)
        dimension x(100000),whts(100000)
        dimension u(100,100000),v(100,100000)
        dimension xs(100000),ws(100000)
c
c       ... orders for which logquad_aux.f actually carries a
c       logquad2xN block.  hqsuppauxquad has no else-branch, so an
c       unsupported order would silently return an empty rule.
c
        dimension navail(14)
        data navail /1,2,3,4,5,6,8,10,12,14,16,18,20,24/
c
        call prini(0,0)
c
 1000   continue
        read(*,*,end=9000) n
c
        iok = 0
        do i=1,14
        if( n .eq. navail(i) ) iok = 1
        enddo
        if( iok .eq. 0 ) then
        write(0,*) 'gentab: no logquad2xN table for order', n
        stop 1
        endif
c
c       ... auxiliary target nodes.  itype=1 returns before the
c       trailing [0,1] remap, so x/whts come back on [-1,1].
c
        itype = 1
        nquad = n
        call legeexps_log_lr(itype,n,nquad,m,x,u,v,whts)
c
        write(*,'(a,1x,i4,1x,i6)') 'ORDER', n, m
        write(*,'(a)') 'AUX'
        do i=1,m
        write(*,'(2(1x,e42.34))') x(i), whts(i)
        enddo
c
c       ... per-target source rules.  hqsuppauxquad returns nodes on
c       [0,1]; map to [-1,1] and rescale the weights accordingly.
c       Node order is left as returned (targets in the right half
c       carry mirrored, hence descending, nodes).
c
        do inode=1,m
        npts = 0
        call hqsuppauxquad(n,inode,xs,ws,npts)
        if( npts .le. 0 ) then
        write(0,*) 'gentab: empty rule, order', n, ' inode', inode
        stop 1
        endif
        write(*,'(a,1x,i4,1x,i6)') 'SRC', inode, npts
        do i=1,npts
        write(*,'(2(1x,e42.34))') 2*xs(i)-1, 2*ws(i)
        enddo
        enddo
c
        goto 1000
 9000   continue
c
        stop
        end
