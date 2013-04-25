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

/**
 * @author Guenther Grau
 */
public class DetectionResultColumn {
  private static final int MAX_NEARBY_DISTANCE = 5;
  protected final BoundingBox boundingBox;
  private final Codeword[] codewords;

  public DetectionResultColumn(final BoundingBox boundingBox) {
    this.boundingBox = new BoundingBox(boundingBox);
    codewords = new Codeword[boundingBox.getMaxY() - boundingBox.getMinY() + 1];
  }

  public Codeword getCodewordNearby(int imageRow) {
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

  protected int getCodewordsIndex(int imageRow) {
    return imageRow - boundingBox.getMinY();
  }

  public int getImageRow(int codewordIndex) {
    return boundingBox.getMinY() + codewordIndex;
  }

  public void setCodeword(int imageRow, Codeword codeword) {
    codewords[getCodewordsIndex(imageRow)] = codeword;
  }

  public Codeword getCodeword(int imageRow) {
    return codewords[getCodewordsIndex(imageRow)];
  }

  public BoundingBox getBoundingBox() {
    return boundingBox;
  }

  public Codeword[] getCodewords() {
    return codewords;
  }
}
