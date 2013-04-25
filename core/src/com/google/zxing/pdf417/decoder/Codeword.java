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
public class Codeword {
  protected static final int BARCODE_ROW_UNKNOWN = -1;

  private final int startX;
  private final int endX;
  private final int bucket;
  private final int value;
  private int rowNumber = BARCODE_ROW_UNKNOWN;

  public Codeword(int startX, int endX, int bucket, int value) {
    this.startX = startX;
    this.endX = endX;
    this.bucket = bucket;
    this.value = value;
  }

  public boolean hasValidRowNumber() {
    return isValidRowNumber(rowNumber);
  }

  public boolean isValidRowNumber(int rowNumber) {
    return BARCODE_ROW_UNKNOWN != rowNumber && bucket == (rowNumber % 3) * 3;
  }

  public void setRowNumberAsRowIndicatorColumn() {
    rowNumber = (value / 30) * 3 + bucket / 3;
  }

  public int getWidth() {
    return endX - startX;
  }

  public int getStartX() {
    return startX;
  }

  public int getEndX() {
    return endX;
  }

  public int getBucket() {
    return bucket;
  }

  public int getValue() {
    return value;
  }

  public int getRowNumber() {
    return rowNumber;
  }

  public void setRowNumber(int rowNumber) {
    this.rowNumber = rowNumber;
  }
}
