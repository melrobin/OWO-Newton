#define F77_FUNC(x) x##_

int add_(int x, int y)
{
  return(x+y);
}
int main()
{
  add_(2,3);
  F77_FUNC(add)(2,3);
  return 0;
}

