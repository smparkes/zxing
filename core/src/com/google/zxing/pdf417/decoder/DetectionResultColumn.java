/*
 * Copyright 2013 ZXing authors
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

import java.util.Formatter;

/**
 * @author Guenther Grau
 */
class DetectionResultColumn {
  
  private static final int MAX_NEARBY_DISTANCE = 5;
  
  private final BoundingBox boundingBox;
  private final Codeword[] codewords;

  DetectionResultColumn(BoundingBox boundingBox) {
    this.boundingBox = new BoundingBox(boundingBox);
    codewords = new Codeword[boundingBox.getMaxY() - boundingBox.getMinY() + 1];
  }

  final Codeword getCodewordNearby(int imageRow) {
    Codeword codeword = getCodeword(imageRow);
    if (codeword != null) {
      return codeword;
    }
    for (int i = 1; i < MAX_NEARBY_DISTANCE; i++) {
      int nearImageRow = getCodewordsIndex(imageRow) - i;
      if (nearImageRow >= 0) {
        codeword = codewords[nearImageRow];
        if (codeword != null) {
          return codeword;
        }
      }
      nearImageRow = getCodewordsIndex(imageRow) + i;
      if (nearImageRow < codewords.length) {
        codeword = codewords[nearImageRow];
        if (codeword != null) {
          return codeword;
        }
      }
    }
    return null;
  }

  final int getCodewordsIndex(int imageRow) {
    return imageRow - boundingBox.getMinY();
  }

  final int getImageRow(int codewordIndex) {
    return boundingBox.getMinY() + codewordIndex;
  }

  final void setCodeword(int imageRow, Codeword codeword) {
    codewords[getCodewordsIndex(imageRow)] = codeword;
  }

  final Codeword getCodeword(int imageRow) {
    return codewords[getCodewordsIndex(imageRow)];
  }

  final BoundingBox getBoundingBox() {
    return boundingBox;
  }

  final Codeword[] getCodewords() {
    return codewords;
  }

  public String getLogString() {
    Formatter formatter = new Formatter();
    int row = 0;
    for (Codeword codeword : codewords) {
      if (codeword == null) {
        formatter.format("%3d:    |   \n", row++);
        continue;
      }
      formatter.format("%3d: %3d|%3d\n", row++, codeword.getRowNumber(), codeword.getValue());
    }
    String result = formatter.toString();
    formatter.close();
    return result;
  }

}
