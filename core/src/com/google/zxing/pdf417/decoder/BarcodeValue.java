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

import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;

/**
 * @author Guenther Grau
 */
public class BarcodeValue {
  Map<Integer,Integer> values = new HashMap<Integer,Integer>();

  public void setValue(int value) {
    Integer confidence = values.get(value);
    if (confidence == null) {
      confidence = 0;
    }
    confidence = confidence + 1;
    values.put(value, confidence);
  }

  public Integer getValue() {
    int maxConfidence = -1;
    Integer result = null;
    boolean ambigous = false;
    for (Entry<Integer,Integer> entry : values.entrySet()) {
      if (entry.getValue() > maxConfidence) {
        maxConfidence = entry.getValue();
        result = entry.getKey();
        ambigous = false;
      } else if (entry.getValue() > maxConfidence) {
        ambigous = true;
      }
    }
    if (ambigous) {
      return null;
    }
    return result;
  }

  public Integer getConfidence(int value) {
    return values.get(value);
  }
}
