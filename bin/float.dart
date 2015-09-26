// Copyright (c) 2015, the Dart project authors.  Please see the AUTHORS file
// for details. All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

library float;

import 'dart:typed_data';  // Only used for converting a double to its bits.

class DoubleProperties {
  final int bitRepresentation;
  final double value;

  static const int exponentMask = 0x7FF0000000000000;
  static const int significandMask = 0x000FFFFFFFFFFFFF;
  static const int signMask = 0x8000000000000000;
  static const int hiddenBit = 0x0010000000000000;
  static const int physicalSignificandSize = 52;
  static const int significandSize = 53;

  static const int exponentBias = 0x3FF + physicalSignificandSize;
  static const int denormalExponent = -exponentBias + 1;

  DoubleProperties(double value)
      : bitRepresentation = _doubleToBitRepresentation(value),
        this.value = value;

  DoubleProperties.fromBits(int bits)
      : bitRepresentation = bits,
        value = _bitRepresentationToDouble(bits);

  static int _doubleToBitRepresentation(double d) {
    // TODO(floitsch): this could be much cheaper if there was a native call.
    Float64List floatList = new Float64List(1);
    // TODO(floitsch): maybe we should care little-endian/big-endian here.
    Int64List int64List = new Int64List.view(floatList.buffer);
    floatList[0] = d;
    return int64List[0];
  }

  static double _bitRepresentationToDouble(int bits) {
    // TODO(floitsch): this could be much cheaper if there was a native call.
    Float64List floatList = new Float64List(1);
    // TODO(floitsch): maybe we should care little-endian/big-endian here.
    Int64List int64List = new Int64List.view(floatList.buffer);
    int64List[0] = bits;
    return floatList[0];
  }

  bool get isDenormal => (bitRepresentation & exponentMask) == 0;

  int get significand {
    int result = bitRepresentation & significandMask;
    if (isDenormal) return result;
    return result | hiddenBit;
  }

  int get exponent {
    if (isDenormal) return denormalExponent;

    int biasedE = (bitRepresentation & exponentMask) >> physicalSignificandSize;
    return biasedE - exponentBias;
  }

  bool get isLowerBoundaryCloser {
    // The boundary is closer if the significand is of the form f == 2^p-1 then
    // the lower boundary is closer.
    // Think of v = 1000e10 and v- = 9999e9.
    // Then the boundary (== (v - v-)/2) is not just at a distance of 1e9 but
    // at a distance of 1e8.
    // The only exception is for the smallest normal: the largest denormal is
    // at the same distance as its successor.
    // Note: denormals have the same exponent as the smallest normals.
    bool physicalSignificandIsZero =
    (bitRepresentation & significandMask) == 0;
    return physicalSignificandIsZero && (exponent != denormalExponent);
  }

  bool get isNegative => (bitRepresentation & signMask) != 0;

  /// Returns the next double. This double must be positive.
  double get nextDouble {
    assert(value != double.INFINITY && !isNegative);
    return new DoubleProperties.fromBits(bitRepresentation + 1).value;
  }

  /// Returns the previous double. This double must be strictly positive.
  double get previousDouble {
    assert(value > 0);
    return new DoubleProperties.fromBits(bitRepresentation - 1).value;
  }

  static const double minPositiveNonZeroDouble = 4e-324;
}
