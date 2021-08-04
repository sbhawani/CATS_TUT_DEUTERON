/* rand example: guess the number */
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

double P(int x, int y, int z) {
  return 2 * z - x + y;
}

void Randum() {
  int x, y, z;

  /* initialize random seed: */
  srand (time(NULL));

  /* generate secret number between 1 and 10: */
  x = rand() % 120 + 0;
  y = rand() % 120 + 0;
  z = rand() % 120 + 0;
  double valP = 0;

  if (x + y + z <= 120) {
    valP = P(x, y, z);

    do {
      printf ("Guess the number (1 to 10): ");
      scanf ("%d", &valP);
      if (iSecret < iGuess) puts ("The secret number is lower");
      else if (iSecret > iGuess) puts ("The secret number is higher");
    } while (iSecret != iGuess);
  }
  puts ("Congratulations!");
  return 0;
}