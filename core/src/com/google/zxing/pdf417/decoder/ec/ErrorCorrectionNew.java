/*
 * Copyright 2012 ZXing authors
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.google.zxing.pdf417.decoder.ec;

import com.google.zxing.ChecksumException;

/**
 * <p>PDF417 error correction implementation.</p>
 *
 * <p>This <a href="http://en.wikipedia.org/wiki/Reed%E2%80%93Solomon_error_correction#Example">example</a>
 * is quite useful in understanding the algorithm.</p>
 *
 * @author Sean Owen
 * @see com.google.zxing.common.reedsolomon.ReedSolomonDecoder
 */
public final class ErrorCorrectionNew {

  private final ModulusGF field;

  public ErrorCorrectionNew() {
    this.field = ModulusGF.PDF417_GF;
  }

  /**
   * @return number of errors
   */
  public int decode_old(int[] received, int numECCodewords, int[] erasures) throws ChecksumException {

    ModulusPoly poly = new ModulusPoly(field, received);
    int[] S = new int[numECCodewords];
    boolean error = false;
    for (int i = numECCodewords; i > 0; i--) {
      int eval = poly.evaluateAt(field.exp(i));
      S[numECCodewords - i] = eval;
      if (eval != 0) {
        error = true;
      }
    }

    if (!error) {
      return 0;
    }

    ModulusPoly knownErrors = field.getOne();
    for (int erasure : erasures) {
      int b = field.exp(received.length - 1 - erasure);
      // Add (1 - bx) term:
      ModulusPoly term = new ModulusPoly(field, new int[] {field.subtract(0, b), 1});
      knownErrors = knownErrors.multiply(term);
    }

    ModulusPoly syndrome = new ModulusPoly(field, S);
    //syndrome = syndrome.multiply(knownErrors);

    ModulusPoly[] sigmaOmega = runEuclideanAlgorithm(field.buildMonomial(numECCodewords, 1), syndrome, numECCodewords);
    //ModulusPoly[] sigmaOmega = runExtendedEuclideanAlgorithm(syndrome, numECCodewords, erasures);
    ModulusPoly sigma = sigmaOmega[0];
    ModulusPoly omega = sigmaOmega[1];

    //sigma = sigma.multiply(knownErrors);

    int[] errorLocations = findErrorLocations(sigma);
    int[] errorMagnitudes = findErrorMagnitudes(omega, sigma, errorLocations);

    for (int i = 0; i < errorLocations.length; i++) {
      int position = received.length - 1 - field.log(errorLocations[i]);
      if (position < 0) {
        throw ChecksumException.getChecksumInstance();
      }
      received[position] = field.subtract(received[position], errorMagnitudes[i]);
    }
    return errorLocations.length;
  }

  private static final int NN = 1024;
  private static final int PRIM = 1;
  private static final int GPRIME = 929;
  private static final int KK = 32;
  private static final int A0 = 928;
  private static final int Ldec = 1;

  int[] Alpha_to = new int[1024];
  int[] Index_of = new int[1024];

  static boolean rs_init = false;

  // initialize table of 3**i for syndrome calculation
  void powers_init() {
    int ii;
    int power_of_3;
    int debug;

    debug = 0;
    power_of_3 = 1;
    Index_of[1] = GPRIME - 1;

    for (ii = 0; ii < GPRIME - 1; ii += 1) {
      Alpha_to[ii] = power_of_3;
      if (power_of_3 < GPRIME) {
        if (ii != GPRIME - 1)
          Index_of[power_of_3] = ii;
      } else {
        System.err.println("Internal error: powers of 3 calculation\n");
      }
      /*
            if (debug) {
                printf("pow = %d  ii = %d \n", power_of_3, ii);
            }
            if (debug) {
                printf("log = %d \n", ii);
            }
            */
      power_of_3 = (power_of_3 * 3) % GPRIME;
    }
    Index_of[0] = GPRIME - 1;
    Alpha_to[GPRIME - 1] = 1;
    Index_of[GPRIME] = A0;
  }

  private int modbase(int x) {
    return ((x) % (GPRIME - 1));
  }

  public int decode(int[] received, int numECCodewords, int[] erasures) throws ChecksumException {
    int[] data = new int[NN];
    System.arraycopy(received, 0, data, 0, received.length);
    int[] eras_pos = new int[NN];
    System.arraycopy(erasures, 0, eras_pos, 0, erasures.length);
    int errorCount = eras_dec_rs(data, eras_pos, erasures.length, received.length, numECCodewords);
    System.arraycopy(data, 0, received, 0, received.length);
    return errorCount;
  }

  /*
   * Performs ERRORS+ERASURES decoding of RS codes. If decoding is successful,
   * writes the codeword into data[] itself. Otherwise data[] is unaltered.
   *
   * Return number of symbols corrected, or -1 if codeword is illegal
   * or uncorrectable. If eras_pos is non-null, the detected error locations
   * are written back. NOTE! This array must be at least NN-KK elements long.
   * 
   * First "no_eras" erasures are declared by the calling program. Then, the
   * maximum # of errors correctable is t_after_eras = floor((NN-KK-no_eras)/2).
   * If the number of channel errors is not greater than "t_after_eras" the
   * transmitted codeword will be recovered. Details of algorithm can be found
   * in R. Blahut's "Theory ... of Error-Correcting Codes".
   *
   * Warning: the eras_pos[] array must not contain duplicate entries; decoder failure
   * will result. The decoder *could* check for this condition, but it would involve
   * extra time on every decoding operation.
   */

  //  int eras_dec_rs(int data[NN], int eras_pos[NN - KK], int no_eras,
  //                  int data_len, int synd_len)
  int eras_dec_rs(int[] data, int[] eras_pos, int no_eras, int data_len, int synd_len) {
    int deg_lambda, el, deg_omega;
    int i, j, r, k;
    int u, q, tmp, num1, num2, den, discr_r;
    int[] lambda = new int[2048 + 1];
    int[] s = new int[2048 + 1]; // Err+Eras Locator poly and syndrome poly
    int[] b = new int[2048 + 1];
    int[] t = new int[2048 + 1];
    int[] omega = new int[2048 + 1];
    int[] root = new int[2048];
    int[] reg = new int[2048 + 1];
    int[] loc = new int[2048];
    boolean syn_error;
    int count;
    int ci;
    int error_val;
    int fix_loc;
    boolean debug;

    if (!rs_init)
      powers_init();

    /* form the syndromes; i.e. evaluate data(x) at roots of g(x)
       namely @**(1+i)*PRIM, i = 0, ... , (NN-KK-1) */
    for (i = 1; i <= synd_len; i++) {
      s[i] = 0;//data[data_len];
    }

    for (j = 1; j <= data_len; j++) {

      if (data[data_len - j] == 0)
        continue;

      tmp = Index_of[data[data_len - j]];

      /*  s[i] ^= Alpha_to[modbase(tmp + (1+i-1)*j)]; */
      for (i = 1; i <= synd_len; i++) {
        s[i] = (s[i] + Alpha_to[modbase(tmp + (i) * j)]) % GPRIME;
      }
    }

    /* Convert syndromes to index form, checking for nonzero condition */
    syn_error = false;
    for (i = 1; i <= synd_len; i++) {
      syn_error |= s[i] != 0;
      //      if (debug) {
      //          System.err.println("Raw syndrome = %d i = %d \n", s[i], i);
      //      }
      s[i] = Index_of[s[i]];
    }

    if (!syn_error) {
      /* if syndrome is zero, data[] is a codeword and there are no
       * errors to correct. So return data[] unmodified
       */
      //      printf("No errors \n");
      return 0;
    }

    for (ci = synd_len - 1; ci >= 0; ci--)
      lambda[ci + 1] = 0;

    lambda[0] = 1;

    if (no_eras > 0) {
      /* Init lambda to be the erasure locator polynomial */
      lambda[1] = Alpha_to[modbase(PRIM * eras_pos[0])];
      for (i = 1; i < no_eras; i++) {
        u = modbase(PRIM * eras_pos[i]);
        for (j = i + 1; j > 0; j--) {
          tmp = Index_of[lambda[j - 1]];
          if (tmp != A0)
            lambda[j] = (lambda[j] + Alpha_to[modbase(u + tmp)]) % GPRIME;
        }
      }
      //  #if DEBUG >= 1
      //      /* Test code that verifies the erasure locator polynomial just constructed
      //         Needed only for decoder debugging. */
      //
      //      /* find roots of the erasure location polynomial */
      //      for (i = 1; i <= no_eras; i++)
      //          reg[i] = Index_of[lambda[i]];
      //      count = 0;      // usg NN = GPRIME-1?
      //
      //      for (i = 1, k = data_len - Ldec; i <= data_len + synd_len;
      //           i++, k = modbase(data_len + k - Ldec)) {
      //          q = 1;
      //          for (j = 1; j <= no_eras; j++)
      //          if (reg[j] != A0) {
      //              reg[j] = modbase(reg[j] + j);
      //              q = (q + Alpha_to[reg[j]]) % GPRIME;
      //          }
      //          if (q != 0)
      //          continue;
      //          /* store root and error location number indices */
      //          root[count] = i;
      //          loc[count] = k;
      //          count++;
      //      }
      //      if (count != no_eras) {
      //          // printf("\n lambda(x) is WRONG\n");
      //          // count = -1;
      //          //  goto finish;
      //      }
      //  #if DEBUG >= 2
      //      printf
      //          ("\n Erasure positions as determined by roots of Eras Loc Poly:\n");
      //      for (i = 0; i < count; i++)
      //          printf("loc  = %d ", loc[i]);
      //      printf("\n");
      //  #endif
      //  #endif
    }
    for (i = 0; i < synd_len + 1; i++)
      b[i] = Index_of[lambda[i]];

    /*
     * Begin Berlekamp-Massey algorithm to determine error+erasure
     * locator polynomial
     */
    r = no_eras;
    el = no_eras;
    while (++r <= synd_len) { /* r is the step number */
      /* Compute discrepancy at the r-th step in poly-form */
      discr_r = 0;
      for (i = 0; i < r; i++) {
        if ((lambda[i] != 0) && (s[r - i] != A0)) {
          //          if (debug) {
          //              printf("do add Index_of[lambda[]] = %d \n",
          //                 Index_of[lambda[i]]);
          //          }
          if (i % 2 == 1) {
            discr_r = (discr_r + Alpha_to[modbase((Index_of[lambda[i]] + s[r - i]))]) % GPRIME;
          } else {
            discr_r = (discr_r + GPRIME - Alpha_to[modbase((Index_of[lambda[i]] + s[r - i]))]) % GPRIME;
          }
          //          if (debug) {
          //              printf("In loop - discr = %d i = %d r = %d lambda[i] = %d s[r-i] = %d \n",
          //                 discr_r, i, r, lambda[i], s[r - i]);
          //          }
        }
      }
      //      if (debug) {
      //          printf("r = %d Discrepency = %d \n", r, discr_r);
      //      }

      discr_r = Index_of[discr_r]; /* Index form */

      if (discr_r == A0) {
        /* 2 lines below: B(x) <-- x*B(x) */
        //  COPYDOWN(&b[1],b,synd_len);
        //
        //          if (debug) {
        //          printf("Discrepancy = A0\n");
        //          }
        for (ci = synd_len - 1; ci >= 0; ci--)
          b[ci + 1] = b[ci];
        b[0] = A0;
      } else {
        /* 7 lines below: T(x) <-- lambda(x) - discr_r*x*b(x) */
        /*  the T(x) will become the next lambda */

        t[0] = lambda[0];
        for (i = 0; i < synd_len; i++) {
          //          if (debug) {
          //              printf("i = %d b[i] = %d \n", i, b[i]);
          //          }
          if (b[i] != A0) {

            //  t[i+1] =  (lambda[i+1] + GPRIME -
            //              Alpha_to[modbase(discr_r + GPRIME - 1 -  b[i])]) % GPRIME;
            t[i + 1] = (lambda[i + 1] + Alpha_to[modbase(discr_r + b[i])]) % GPRIME;

            //              if (debug) {
            //              printf("New t[i+1] = %d lambda[i+1] = %d b[i] = %d i = %d discr_r = %d\n",
            //                     t[i + 1], lambda[i + 1], b[i], i, discr_r);
            //              }
          } else {
            t[i + 1] = lambda[i + 1];
          }
          //          if (debug) {
          //              printf("i = %d t[i+1] = %d lambda[i+1] = %d \n", i,
          //                 t[i + 1], lambda[i + 1]);
          //          }
        }
        el = 0;
        if (2 * el <= r + no_eras - 1) {
          //          if (debug) {
          //              printf("Reached the el stuff, inv  el = %d r = %d \n", el, r);
          //          }
          el = r + no_eras - el;
          /*
           * 2 lines below: B(x) <-- inv(discr_r) *
           * lambda(x)
           */
          for (i = 0; i <= synd_len; i++) {

            if (lambda[i] == 0) {
              b[i] = A0;
            } else {
              b[i] = modbase(Index_of[lambda[i]] - discr_r + GPRIME - 1);
              //              if (debug) {
              //                  printf("Inverting le  b[i] = %d i = %d \n", b[i], i);
              //              }
            }
          }

        } else {
          //          if (debug) {
          //              printf("Reached the el stuff, x mul,   el = %d r = %d \n", el, r);
          //          }
          /* 2 lines below: B(x) <-- x*B(x) */
          //      COPYDOWN(&b[1],b,synd_len);
          for (ci = synd_len - 1; ci >= 0; ci--)
            b[ci + 1] = b[ci];
          b[0] = A0;
        }
        //      COPY(lambda,t,synd_len+1);

        for (ci = synd_len + 1 - 1; ci >= 0; ci--) {
          lambda[ci] = t[ci];
          //          if (debug) {
          //              printf("ci = %d Lambda = %d \n", ci, t[ci]);
          //          }
        }
      }
    }

    /* Convert lambda to index form and compute deg(lambda(x)) */
    deg_lambda = 0;
    for (i = 0; i < synd_len + 1; i++) {

      lambda[i] = Index_of[lambda[i]];

      if (lambda[i] != A0)
        deg_lambda = i;

      //      if (debug) {
      //          printf("Lambda in index form = %d \n", lambda[i]);
      //      }

    }

    //      if (debug) {
    //      printf("Determination of deg_lambda = %d \n", deg_lambda);
    //      }

    /*
     * Find roots of the error+erasure locator polynomial by Chien
     * Search
     */

    for (ci = synd_len - 1; ci >= 0; ci--)
      reg[ci + 1] = lambda[ci + 1];

    count = 0; /* Number of roots of lambda(x) */
    for (i = 1, k = data_len - 1; i <= GPRIME; i++) {
      q = 1;
      //      if (debug) {
      //          printf(" Reg[j] = %d q = %d i = %d \n", reg[j], q, i);
      //      }
      for (j = deg_lambda; j > 0; j--) {

        if (reg[j] != A0) {
          //          if (debug) {
          //              printf("loop Reg[j] pre = %d \n", reg[j]);
          //          }
          reg[j] = modbase(reg[j] + j);
          //      q = modbase( q +  Alpha_to[reg[j]]);
          if (deg_lambda != 1) {
            if (j % 2 == 0) {
              q = (q + Alpha_to[reg[j]]) % GPRIME;
            } else {
              q = (q + GPRIME - Alpha_to[reg[j]]) % GPRIME;
            }
          } else {
            q = Alpha_to[reg[j]] % GPRIME;
            if (q == 1)
              --q;
          }
          //          if (debug) {
          //              printf("loop Reg[j] = %d q = %d i = %d j = %d %d = k\n",
          //                 reg[j], q, i, j, k);
          //          }
        }
      }

      if (q == 0) {
        /* store root (index-form) and error location number */
        root[count] = i;

        loc[count] = GPRIME - 1 - i;
        if (count < synd_len) {
          count += 1;
        } else {
          System.out.println("Error : Error count too big = %d \n" + count);
        }

        //          if (debug) {
        //          printf("root  = %d loc = %d \n", i, k);
        //          }

      }
      if (k == 0) {
        k = data_len - 1;
      } else {
        k -= 1;
      }

      /* If we've already found max possible roots,
       * abort the search to save time
       */

      if (count == deg_lambda)
        break;

    }

    if (deg_lambda != count) {
      /*
       * deg(lambda) unequal to number of roots => uncorrectable
       * error detected
       */

      System.err.println("Uncorrectable error: root count = %d deg lambda = %d \n" + count + "," + deg_lambda);
      return -1;
    }

    /*
     * Compute err+eras evaluator poly omega(x) = s(x)*lambda(x) (modulo
     * x**(synd_len)). in index form. Also find deg(omega).
     */
    deg_omega = 0;
    for (i = 0; i < synd_len; i++) {
      tmp = 0;
      j = (deg_lambda < i) ? deg_lambda : i;
      //      if (debug) {
      //          printf("j = %d deg_lambda = %d lambda[j] = %d \n",
      //             j, deg_lambda, lambda[j]);
      //      }
      for (; j >= 0; j--) {
        if ((s[i + 1 - j] != A0) && (lambda[j] != A0)) {
          if (j % 2 == 1) {
            tmp = (tmp + GPRIME - Alpha_to[modbase(s[i + 1 - j] + lambda[j])]) % GPRIME;
          } else {

            tmp = (tmp + Alpha_to[modbase(s[i + 1 - j] + lambda[j])]) % GPRIME;
          }
          //          if (debug) {
          //              printf("In tmp loop  tmp = %d i = %d j = %d s[i+1-j] = %d lambda[j] = %d \n",
          //                 tmp, i, j, s[i + 1 - j], lambda[j]);
          //          }
        }
      }

      if (tmp != 0)
        deg_omega = i;
      omega[i] = Index_of[tmp];
      //      if (debug) {
      //          printf("Omega [i] = %d i = %d \n", omega[i], i);
      //      }

    }
    omega[synd_len] = A0;
    //      if (debug) {
    //      printf("Degree of omega = %d \n", deg_omega);
    //      }

    /*
     * Compute error values in poly-form. num1 = omega(inv(X(l))), num2 =
     * inv(X(l))**(B0-1) and den = lambda_pr(inv(X(l))) all in poly-form
     */
    for (j = count - 1; j >= 0; j--) {
      num1 = 0;
      for (i = deg_omega; i >= 0; i--) {
        if (omega[i] != A0) {
          //    num1  = ( num1 + Alpha_to[modbase(omega[i] + (i * root[j])]) % GPRIME;
          num1 = (num1 + Alpha_to[modbase(omega[i] + ((i + 1) * root[j]))]) % GPRIME;
          //          if (debug) {
          //              printf("Num1 = %d i = %d omega[i] = %d root[j] = %d \n",
          //                 num1, i, omega[i], root[j]);
          //          }
        }
      }
      //  num2 = Alpha_to[modbase(root[j] * (1 - 1) + data_len)];

      num2 = 1;
      den = 0;

      // denominator if product of all (1 - Bj Bk) for k != j
      // if count = 1, then den = 1

      den = 1;
      for (k = 0; k < count; k += 1) {
        if (k != j) {
          tmp = (1 + GPRIME - Alpha_to[modbase(GPRIME - 1 - root[k] + root[j])]) % GPRIME;
          den = Alpha_to[modbase(Index_of[den] + Index_of[tmp])];
        }
      }

      //      if (debug) {
      //          printf("den = %d \n", den);
      //      }

      if (den == 0) {
        //  #if DEBUG >= 1
        //          printf("\n ERROR: denominator = 0\n");
        //  #endif
        /* Convert to dual- basis */
        return -1;
      }

      //      if (debug) {
      //          printf("Index num1 = %d Index num2 = %d Index of den = %d \n",
      //             Index_of[num1], Index_of[num2], Index_of[den]);
      //      }

      error_val = Alpha_to[modbase(Index_of[num1] + Index_of[num2] + GPRIME - 1 - Index_of[den])] % GPRIME;

      //      if (debug) {
      //          printf("error_val = %d \n", error_val);
      //      }

      /* Apply error to data */
      if (num1 != 0) {
        if (loc[j] < data_len + 1) {
          fix_loc = data_len - loc[j];
          //fix_loc = loc[j];
          //          if (debug) {
          //              printf("Fix loc = %d \n", fix_loc);
          //          }
          if (fix_loc < data_len + 1) {
            data[fix_loc] = (data[fix_loc] + GPRIME - error_val) % GPRIME;
          }
        }
      }
    }
    //      if (debug) {
    //      printf("At FINISH \n");
    //      }

    if (eras_pos != null) {
      for (int zz = 0; zz < count; zz++) {
        eras_pos[zz] = loc[zz];
      }
    }
    return count;
  }

  private ModulusPoly[] runExtendedEuclideanAlgorithm(ModulusPoly v, int twoT, int[] erasures) throws ChecksumException {
    int delta = -1;
    int[] uCoefficients = new int[twoT + 1];
    uCoefficients[0] = 1;
    ModulusPoly u = new ModulusPoly(field, uCoefficients);
    ModulusPoly x = field.getOne();
    ModulusPoly w = field.getZero();
    ModulusPoly psi = new ModulusPoly(field, erasures);

    ModulusPoly zett = new ModulusPoly(field, new int[] {1, 0});

    for (int i = 0; i < twoT; i++) {
      boolean first = psi.evaluateAt(0) != 0;
      boolean swap = !first && (v.evaluateAt(twoT - 1) != 0) && delta < 0;
      int gamma = first ? psi.evaluateAt(0) : u.evaluateAt(twoT);
      int xi = first ? 1 : v.evaluateAt(twoT - 1);
      if (swap) {
        delta = -delta - 1;
      } else if (!first) {
        delta--;
      }
      psi = psi.divide(zett)[0];
      v = v.multiply(gamma).multiply(zett).subtract(first ? v.multiply(xi) : u.multiply(xi));
      x = x.multiply(gamma).multiply(zett).subtract(first ? x.multiply(xi) : w.multiply(xi));
      if (swap) {
        u = v.multiply(zett);
        w = x.multiply(zett);
      }
    }
    if (delta >= 0 || psi.evaluateAt(0) != 0) {
      throw ChecksumException.getChecksumInstance();
    }

    return new ModulusPoly[] {x, v};
  }

  private ModulusPoly[] runEuclideanAlgorithm(ModulusPoly a, ModulusPoly b, int R) throws ChecksumException {
    // Assume a's degree is >= b's
    if (a.getDegree() < b.getDegree()) {
      ModulusPoly temp = a;
      a = b;
      b = temp;
    }

    ModulusPoly rLast = a;
    ModulusPoly r = b;
    ModulusPoly tLast = field.getZero();
    ModulusPoly t = field.getOne();

    // Run Euclidean algorithm until r's degree is less than R/2
    while (r.getDegree() >= R / 2) {
      ModulusPoly rLastLast = rLast;
      ModulusPoly tLastLast = tLast;
      rLast = r;
      tLast = t;

      // Divide rLastLast by rLast, with quotient in q and remainder in r
      if (rLast.isZero()) {
        // Oops, Euclidean algorithm already terminated?
        throw ChecksumException.getChecksumInstance();
      }
      r = rLastLast;
      ModulusPoly q = field.getZero();
      int denominatorLeadingTerm = rLast.getCoefficient(rLast.getDegree());
      int dltInverse = field.inverse(denominatorLeadingTerm);
      while (r.getDegree() >= rLast.getDegree() && !r.isZero()) {
        int degreeDiff = r.getDegree() - rLast.getDegree();
        int scale = field.multiply(r.getCoefficient(r.getDegree()), dltInverse);
        q = q.add(field.buildMonomial(degreeDiff, scale));
        r = r.subtract(rLast.multiplyByMonomial(degreeDiff, scale));
      }

      t = q.multiply(tLast).subtract(tLastLast).negative();
    }

    int sigmaTildeAtZero = t.getCoefficient(0);
    if (sigmaTildeAtZero == 0) {
      throw ChecksumException.getChecksumInstance();
    }

    int inverse = field.inverse(sigmaTildeAtZero);
    ModulusPoly sigma = t.multiply(inverse);
    ModulusPoly omega = r.multiply(inverse);
    return new ModulusPoly[] {sigma, omega};
  }

  private int[] findErrorLocations(ModulusPoly errorLocator) throws ChecksumException {
    // This is a direct application of Chien's search
    int numErrors = errorLocator.getDegree();
    int[] result = new int[numErrors];
    int e = 0;
    for (int i = 1; i < field.getSize() && e < numErrors; i++) {
      if (errorLocator.evaluateAt(i) == 0) {
        result[e] = field.inverse(i);
        e++;
      }
    }
    if (e != numErrors) {
      throw ChecksumException.getChecksumInstance();
    }
    return result;
  }

  private int[] findErrorMagnitudes(ModulusPoly errorEvaluator, ModulusPoly errorLocator, int[] errorLocations) {
    int errorLocatorDegree = errorLocator.getDegree();
    int[] formalDerivativeCoefficients = new int[errorLocatorDegree];
    for (int i = 1; i <= errorLocatorDegree; i++) {
      formalDerivativeCoefficients[errorLocatorDegree - i] = field.multiply(i, errorLocator.getCoefficient(i));
    }
    ModulusPoly formalDerivative = new ModulusPoly(field, formalDerivativeCoefficients);

    // This is directly applying Forney's Formula
    int s = errorLocations.length;
    int[] result = new int[s];
    for (int i = 0; i < s; i++) {
      int xiInverse = field.inverse(errorLocations[i]);
      int numerator = field.subtract(0, errorEvaluator.evaluateAt(xiInverse));
      int denominator = field.inverse(formalDerivative.evaluateAt(xiInverse));
      result[i] = field.multiply(numerator, denominator);
    }
    return result;
  }
}
