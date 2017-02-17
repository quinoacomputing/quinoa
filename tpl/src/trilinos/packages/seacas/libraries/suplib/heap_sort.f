      subroutine heapsort(a,n)
      
      dimension a(0:*)
      integer n
      
      integer start, bottom
      real temp
 
      do start = (n-2)/2, 0, -1
        call hs_siftdown(a, start, n);
      end do
 
      do bottom = n-1, 1, -1
        temp = a(0)
        a(0) = a(bottom)
        a(bottom) = temp;
        call hs_siftdown(a, 0, bottom)
      end do
      end
 
      subroutine hs_siftdown(a, start, bottom)
 
      real a(0:*)
      integer start, bottom
      integer child, root
      real temp
 
      root = start
      do while(root*2 + 1 < bottom)
        child = root * 2 + 1
 
        if ((child + 1 < bottom) .and. (a(child) < a(child+1))) then
          child = child + 1
        end if
 
        if (a(root) < a(child)) then
          temp = a(child)
          a(child) = a (root)
          a(root) = temp
          root = child
        else
          return
        end if  
      end do    
      return
      end


