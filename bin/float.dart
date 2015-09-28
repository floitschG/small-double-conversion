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


  static const double minPositiveNonZero = 4e-324;

  /// The greatest decimal exponent that can be encoded in a double without
  /// overflowing to infinity: 1e308.
  static const int maxDecimalExponent = 308;

  /// 1/10^maxDecimalExponent.
  static const double invMaxDecimalPower = 1e-308;

  /// A decimal of size 15 fits into a double without loss of precision.
  ///
  /// The significant of a double is 53 bits (when counting the hidden bit).
  ///
  ///     2^53 = 9007199254740992.
  ///
  /// Any integer with at most 15 decimal digits will hence fit into a double
  /// without loss of precision.
  static const int maxExactDecimalDigits = 15;

  // If size is important this table can be replaced by a call to
  // `math.pow(2.0, i)`.
  static const List<double> exactPowersOfTen = const <double>[
    1.0,  // 10^0
    10.0,
    100.0,
    1000.0,
    10000.0,
    100000.0,
    1000000.0,
    10000000.0,
    100000000.0,
    1000000000.0,
    // 10^10 = 0x2540be400 = 0x9502f9 * 2^10
    10000000000.0,  // 10^10
    100000000000.0,
    1000000000000.0,
    10000000000000.0,
    100000000000000.0,
    1000000000000000.0,
    10000000000000000.0,
    100000000000000000.0,
    1000000000000000000.0,
    10000000000000000000.0,
    100000000000000000000.0,  // 10^20
    1000000000000000000000.0,
    // 10^22 = 0x21e19e0c9bab2400000 = 0x878678326eac9 * 2^22
    10000000000000000000000.0
  ];

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
  double get next {
    assert(value != double.INFINITY && !isNegative);
    return new DoubleProperties.fromBits(bitRepresentation + 1).value;
  }

  /// Returns the previous double. This double must be strictly positive.
  double get previous {
    assert(value > 0);
    return new DoubleProperties.fromBits(bitRepresentation - 1).value;
  }
}

class SingleProperties {
  final int bitRepresentation;
  final double value;

  static const int exponentMask = 0x7F800000;
  static const int significandMask = 0x007FFFFF;
  static const int signMask = 0x80000000;
  static const int hiddenBit = 0x00800000;
  static const int physicalSignificandSize = 23;
  static const int significandSize = 24;

  static const int exponentBias = 0x7F + physicalSignificandSize;
  static const int denormalExponent = -exponentBias + 1;

  /// The smallest positive non-zero single.
  static const double minPositiveNonZero = 1e-45;

  /// The greatest decimal exponent that can be encoded in a single without
  /// overflowing to infinity: 1e38.
  static const int maxDecimalExponent = 38;

  /// 1/10^maxDecimalExponent.
  static const double invMaxDecimalPower = 1e-38;

  /// A decimal of size 7 fits into a single without loss of precision.
  ///
  /// The significant of a single is 24 bits (when counting the hidden bit).
  ///
  ///     2^24 = 16777216.
  ///
  /// Any integer with at most 7 decimal digits will hence fit into a single
  /// without loss of precision.
  static const int maxExactDecimalDigits = 7;

  // If size is important this table can be replaced by a call to
  // `math.pow(2.0, i)`.
  static const List<double> exactPowersOfTen = const <double>[
    1.0,  // 10^0
    10.0,
    100.0,
    1000.0,
    10000.0,
    100000.0,
    1000000.0,
    10000000.0,
    100000000.0,
    1000000000.0,
    // 10^10 = 0x2540be400 = 0x9502f9 * 2^10
    10000000000.0,  // 10^10
  ];


  SingleProperties(double value)
      : bitRepresentation = _singleToBitRepresentation(value),
        this.value = value;

  SingleProperties.fromBits(int bits)
      : bitRepresentation = bits,
        value = _bitRepresentationToSingle(bits);

  static int _singleToBitRepresentation(double d) {
    // TODO(floitsch): this could be much cheaper if there was a native call.
    Float32List floatList = new Float32List(1);
    // TODO(floitsch): maybe we should care little-endian/big-endian here.
    Int32List int32List = new Int32List.view(floatList.buffer);
    floatList[0] = d;
    return int32List[0];
  }

  static double _bitRepresentationToSingle(int bits) {
    // TODO(floitsch): this could be much cheaper if there was a native call.
    Float32List floatList = new Float32List(1);
    // TODO(floitsch): maybe we should care little-endian/big-endian here.
    Int32List int32List = new Int32List.view(floatList.buffer);
    int32List[0] = bits;
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
  double get next {
    assert(value != double.INFINITY && !isNegative);
    return new DoubleProperties.fromBits(bitRepresentation + 1).value;
  }

  /// Returns the previous double. This double must be strictly positive.
  double get previous {
    assert(value > 0);
    return new DoubleProperties.fromBits(bitRepresentation - 1).value;
  }
}
