/* Copyright (c) 2007 Michele Bini
 * 
 * This program is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "curve25519.h"
#include "base32.h"
#include "tsc_x86.h"

static void
usage(FILE*f, const char*p, int q) {
  if (q) {
    fprintf(f, "Try '%s --help' for more information.\n", p);
    return;
  }
  fprintf(f,
	  "Usage:\n\n"
	  "Generating a public key from your private key:\n"
	  "  %s [OPT]... [FMT] <private key> [FMT]\n\n"
	  "Obtain shared secret from <public key> with your private key:\n"
	  "  %s [OPT]... [FMT] <private key> [FMT] <public key> [FMT]\n\n"
	  "FMT specifies the format used for the keys.\n"
	  "It is one of the options:\n"
	  "  --b32: base32-encoded (default)\n"
	  "  --ibh: inverted-bytes hexadecimal\n"
	  "  --hex: hexadecimal\n\n"
	  "OPT is one of the options:\n"
	  "  --safe: abort program when potentially unsafe keys are seen\n"
	  "  --warn: just warn about it (default)\n"
	  "  --unsafe: do not perform any key validation\n\n",
	  p, p, p);
}

int main(int argc, const char *argv[]) {
  curve25519key_t *k = malloc(sizeof(curve25519key_t));
  int format = 0; /* 0: base32; 1: hex; 2: byte-inverted hex */
  int c = 1;
  int kk = 0; /* number of keys parsed */
  int sf = 1; /* 0: no validation; 1: warn; 2: reject invalid keys */
  while (c<argc) {
    int t = 0;
    const char*a = argv[c];
    c++;
    if (strlen(a) > 8) {
      t = 1;
    } else if (a[0] == '-') {
      if (strcmp(a, "--hex") == 0) {
	format = 1;
      } else if (strcmp(a, "--ibh") == 0) {
	format = 2;
      } else if (strcmp(a, "--b32") == 0) {
	format = 0;
      } else if (strcmp(a, "--safe") == 0) {
	sf = 2;
      } else if (strcmp(a, "--warn") == 0) {
	sf = 1;
      } else if (strcmp(a, "--unsafe") == 0) {
	sf = 0;
      } else if (strcmp(a, "--help") == 0) {
	usage(stdout, argv[0], 0);
	exit(EXIT_SUCCESS);
      }
    } else {
      t = 1;
    }
    if (t == 1) {
      k = realloc(k, sizeof(curve25519key_t)*(kk + 1));
      switch (format) {
      case 0:
	base32_decode(a, k + kk);
	break;

      case 1:
	{
	  int l = strlen(a);
	  unsigned int p;
	  unsigned int b = 0;
	  while (l > 0) {
	    l--; p = a[l];
	    switch (p) {
	    case '0': case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7': case '8': case '9':
	      p = p - '0';
	      break;
	    case 'a': case 'b': case 'c': case 'd': case 'e': case 'f':
	      p = p - ('a' - 10);
	      break;
	    case 'A': case 'B': case 'C': case 'D': case 'E': case 'F':
	      p = p - ('A' - 10);
	      break;
	    default:
	      fprintf(stderr, "Bad character where an hexadecimal key was expected: %c\n", p);
	      exit(EXIT_FAILURE);
	    }
	    curve25519key_setbit(k+kk, b, p&1); p >>= 1; b++;
	    curve25519key_setbit(k+kk, b, p&1); p >>= 1; b++;
	    curve25519key_setbit(k+kk, b, p&1); p >>= 1; b++;
	    curve25519key_setbit(k+kk, b, p&1); p >>= 1; b++;
	  }
	  while (b < C25519BITS) {
	    curve25519key_setbit(k+kk, b, 0); b++;
	  }
	}
	break;

      default:
	{
	  unsigned int x[32];
	  int n = sscanf(a, "%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x",
			 x, x+1, x+2, x+3, x+4, x+5, x+6, x+7,
			 x+8, x+9, x+10, x+11, x+12, x+13, x+14, x+15,
			 x+16, x+17, x+18, x+19, x+20, x+21, x+22, x+23,
			 x+24, x+25, x+26, x+27, x+28, x+29, x+30, x+31
			 );
	  if (n != 32) {
	    fprintf(stderr, "Bad key format.\n");
	    usage(stderr, argv[0], 1);
	    exit(EXIT_FAILURE);
	  }
	  for (n = 0; n < 32; n++) {
	    curve25519key_setbyte(k+kk, n, x[n]);
	  }
	}
      }
      kk++;
    }
  }


  if (sf > 0) {
    int q;
    for (q = 0; q < kk; q++) {
      if (!curve25519key_validate(k + q)) {
	fprintf(stderr, "Input key n. %d may be unsafe!\n", q);
	if (sf > 1) {
	  exit(EXIT_FAILURE);
	}
      }
    }
  }
  
  if (kk == 1) {
    curve25519key_t b = { 9 };
	myInt64 start, end;
	start = start_tsc();
    curve25519(k, k, &b);
	end = stop_tsc( start );
	printf("Runtime: %llu \n", end);
  } else if (kk > 1) {
    int q = kk - 1;
    curve25519key_t *l = k + q;
    while (1) {
      q--;
      curve25519(l, k + q, l);
      if (q <= 0) {
	break;
      }
      if ((sf > 0) && !curve25519key_validate(l)) {
	fprintf(stderr, "Intermediate key may be unsafe!\n");
	if (sf > 1) {
	  exit(EXIT_FAILURE);
	}
      }
    }
  } else {
    usage(stderr, argv[0], 0);
    exit(EXIT_FAILURE);
  }

  k += kk - 1;
  if ((sf > 0) && !curve25519key_validate(k)) {
    fprintf(stderr, "Output key may be unsafe!\n");
    if (sf > 1) {
      exit(EXIT_FAILURE);
    }
  }
  switch (format) {
  case 0:
    {
      char s[(C25519BITS/4)+2];
      base32_encode(s, k); printf("%s\n", s);
    }
    break;

  case 1:
    printf("%08x%08x%08x%08x%08x%08x%08x%08x\n",
	   curve25519key_getuint32(k, 7), curve25519key_getuint32(k, 6), curve25519key_getuint32(k, 5), curve25519key_getuint32(k, 4),
	   curve25519key_getuint32(k, 3), curve25519key_getuint32(k, 2), curve25519key_getuint32(k, 1), curve25519key_getuint32(k, 0));
    break;

  default:
    printf("%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x%02x\n",
	   curve25519key_getbyte(k, 0), curve25519key_getbyte(k, 1), curve25519key_getbyte(k, 2), curve25519key_getbyte(k, 3),
	   curve25519key_getbyte(k, 4), curve25519key_getbyte(k, 5), curve25519key_getbyte(k, 6), curve25519key_getbyte(k, 7),
	   curve25519key_getbyte(k, 8), curve25519key_getbyte(k, 9), curve25519key_getbyte(k, 10), curve25519key_getbyte(k, 11),
	   curve25519key_getbyte(k, 12), curve25519key_getbyte(k, 13), curve25519key_getbyte(k, 14), curve25519key_getbyte(k, 15),
	   curve25519key_getbyte(k, 16), curve25519key_getbyte(k, 17), curve25519key_getbyte(k, 18), curve25519key_getbyte(k, 19),
	   curve25519key_getbyte(k, 20), curve25519key_getbyte(k, 21), curve25519key_getbyte(k, 22), curve25519key_getbyte(k, 23),
	   curve25519key_getbyte(k, 24), curve25519key_getbyte(k, 25), curve25519key_getbyte(k, 26), curve25519key_getbyte(k, 27),
	   curve25519key_getbyte(k, 28), curve25519key_getbyte(k, 29), curve25519key_getbyte(k, 30), curve25519key_getbyte(k, 31));
    break;
  }
  exit(EXIT_SUCCESS);
}
