static unsigned int y = 2463534242U;

unsigned int xorshift (void)
{
   y ^= (y << 13);
   y ^= (y >> 17);
   return y ^= (y << 5);
}
