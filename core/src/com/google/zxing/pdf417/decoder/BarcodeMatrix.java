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
import java.util.HashMap;
import java.util.Map;

/**
 * @author Guenther Grau
 */
final class BarcodeMatrix {

  private final Map<String,BarcodeValue> values = new HashMap<String,BarcodeValue>();
  private int maxRow = -1;
  private int maxColumn = -1;

  private static String getKey(int barcodeRow, int barcodeColumn) {
    return barcodeRow + "," + barcodeColumn;
  }

  void setValue(int row, int column, int value) {
    maxRow = Math.max(maxRow, row);
    maxColumn = Math.max(maxColumn, column);
    String key = getKey(row, column);
    BarcodeValue barcodeValue = values.get(key);
    if (barcodeValue == null) {
      barcodeValue = new BarcodeValue();
      values.put(key, barcodeValue);
    }
    barcodeValue.setValue(value);
  }

  public int[] getValue(int row, int column) {
    BarcodeValue barcodeValue = values.get(getKey(row, column));
    return barcodeValue == null ? null : barcodeValue.getValue();
  }

  @Override
  public String toString() {
    Formatter formatter = new Formatter();
    for (int row = 0; row <= maxRow; row++) {
      formatter.format("Row %2d: ", row);
      for (int column = 0; column <= maxColumn; column++) {
        BarcodeValue barcodeValue = values.get(getKey(row, column));
        if (barcodeValue == null || barcodeValue.getValue() == null) {
          formatter.format("        ", (Object[]) null);
        } else {
          formatter.format("%4d(%2d)", barcodeValue.getValue()[0],
              barcodeValue.getConfidence(barcodeValue.getValue()[0]));
        }
      }
      formatter.format("\n");
    }
    String result = formatter.toString();
    formatter.close();
    return result;
  }
}
