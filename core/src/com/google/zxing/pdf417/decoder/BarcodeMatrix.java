package com.google.zxing.pdf417.decoder;

import java.util.HashMap;
import java.util.Map;

public class BarcodeMatrix {
  Map<String,BarcodeValue> values = new HashMap<String,BarcodeValue>();
  int maxRow = -1;
  int maxColumn = -1;

  private String getKey(int barcodeRow, int barcodeColumn) {
    return barcodeRow + "," + barcodeColumn;
  }

  public void setValue(int row, int column, int value) {
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

  public Integer getValue(int row, int column) {
    BarcodeValue barcodeValue = values.get(getKey(row, column));
    if (barcodeValue == null) {
      return null;
    }
    return barcodeValue.getValue();
  }
}
