/*
 * Copyright 2009 ZXing authors
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

package com.google.zxing.pdf417.decoder;

import com.google.zxing.FormatException;
import com.google.zxing.common.BitMatrix;
import com.google.zxing.pdf417.PDF417Common;

/**
 * <p>
 * This class parses the BitMatrix image into codewords.
 * </p>
 *
 * @author SITA Lab (kevin.osullivan@sita.aero)
 */
public final class BitMatrixParser {

  private final BitMatrix bitMatrix;

  private int leftColumnECData = 0;
  private int rightColumnECData = 0;
  private int eraseCount = 0;
  private int[] erasures;
  private int ecLevel = -1;

  BitMatrixParser(BitMatrix bitMatrix) {
    this.bitMatrix = bitMatrix;
  }

  /**
   * To ensure separability of rows, codewords of consecutive rows belong to
   * different subsets of all possible codewords. This routine scans the
   * symbols in the barcode. When it finds a number of consecutive rows which
   * are the same, it assumes that this is a row of codewords and processes
   * them into a codeword array.
   *
   * @return an array of codewords.
   */
  int[] readCodewords() throws FormatException {
    int height = bitMatrix.getHeight();

    erasures = new int[PDF417Common.MAX_CODEWORDS_IN_BARCODE];

    int[] codewords = new int[PDF417Common.MAX_CODEWORDS_IN_BARCODE];
    int next = 0;
    int rowNumber = 0;
    for (int i = 0; i < height; i++) {
      if (rowNumber >= PDF417Common.MAX_ROWS_IN_BARCODE) {
        // Something is wrong, since we have exceeded
        // the maximum rows in the specification.
        throw FormatException.getFormatInstance();
      }
      // Process Row
      next = processRow(rowNumber, codewords, next);
      rowNumber++;
    }
    erasures = trimArray(erasures, eraseCount);
    return trimArray(codewords, next);
  }

  /**
   * Convert the symbols in the row to codewords.
   * Each PDF417 symbol character consists of four bar elements and four space
   * elements, each of which can be one to six modules wide. The four bar and
   * four space elements shall measure 17 modules in total.
   *
   * @param rowNumber   the current row number of codewords.
   * @param codewords   the codeword array to save codewords into.
   * @param next        the next available index into the codewords array.
   * @return the next available index into the codeword array after processing
   *         this row.
   */
  int processRow(int rowNumber, int[] codewords, int next) throws FormatException {
    int width = bitMatrix.getWidth();
    int columnNumber = 0;
    long symbol = 0;
    for (int i = 0; i < width; i += PDF417Common.MODULES_IN_CODEWORD) {
      for (int mask = PDF417Common.MODULES_IN_CODEWORD - 1; mask >= 0; mask--) {
        if (bitMatrix.get(i + (PDF417Common.MODULES_IN_CODEWORD - 1 - mask), rowNumber)) {
          symbol |= 1L << mask;
        }
      }
      if (columnNumber > 0) {
        int cw = PDF417Common.getCodeword(symbol);
        if (cw < 0 && i < width - PDF417Common.MODULES_IN_CODEWORD) {
          // Skip errors on the Right row indicator column
          if (eraseCount >= erasures.length) {
            throw FormatException.getFormatInstance();
          }
          erasures[eraseCount] = next;
          next++;
          eraseCount++;
        } else {
          if (next >= codewords.length) {
            throw FormatException.getFormatInstance();
          }
          codewords[next++] = cw;
        }
      } else {
        // Left row indicator column
        int cw = PDF417Common.getCodeword(symbol);
        if (ecLevel < 0 && rowNumber % 3 == 1) {
          leftColumnECData = cw;
        }
      }
      symbol = 0;
      columnNumber++;
    }
    if (columnNumber > 1) {
      // Right row indicator column is in codeword[next]
      // Overwrite the last codeword i.e. Right Row Indicator
      --next;
      if (ecLevel < 0 && rowNumber % 3 == 2) {
        rightColumnECData = codewords[next];
        if (rightColumnECData == leftColumnECData && leftColumnECData > 0) {
          ecLevel = (rightColumnECData % 30) / 3;
        }
      }
      codewords[next] = 0;
    }
    return next;
  }

  /**
   * Trim the array to the required size.
   *
   * @param array the array
   * @param size  the size to trim it to
   * @return the new trimmed array
   */
  private static int[] trimArray(int[] array, int size) {
    if (size < 0) {
      throw new IllegalArgumentException();
    }
    int[] a = new int[size];
    System.arraycopy(array, 0, a, 0, size);
    return a;
  }

  /**
   * Returns an array of locations representing the erasures.
   */
  public int[] getErasures() {
    return erasures;
  }

  public int getECLevel() {
    return ecLevel;
  }
}
