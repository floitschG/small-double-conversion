// Copyright (c) 2015, the Dart project authors.  Please see the AUTHORS file
// for details. All rights reserved. Use of this source code is governed by a
// BSD-style license that can be found in the LICENSE file.

import 'dart:math' as math;
import 'char_codes.dart' as codes;
import 'float.dart';

String convertToShortest(double d) {
  if (isSpecial(d)) return handleSpecialValues(d);

  if (d == 0.0) {
    return d.isNegative ? "-0.0" : "0.0";
  }

  DecimalDecomposition decomposition =
      doubleToDecimal(d, DoubleConversionMode.shortest);

  int exponent = decomposition.decimalPointPosition - 1;
  if (lowExponentForDecimal <= exponent && exponent < highExponentForDecimal) {
    int digitsLength = decomposition.decimalDigits.length;
    int decimalPoint = decomposition.decimalPointPosition;
    // To distinguish doubles from integers, we require the shortest
    // representation to always have a least one digit after the decimal point.
    return createDecimalRepresentation(
        decomposition, math.max(1, digitsLength - decimalPoint));
  } else {
    return createExponentialRepresentation(decomposition);
  }
}
String convertToFixed(double d, int requestedDigits) {
  assert(maxFixedDigitsBeforePoint == 60);
  double firstNonFixed = 1e60;

  if (isSpecial(d)) return handleSpecialValues(d);

  if (requestedDigits > maxFixedDigitsAfterPoint) return null;
  if (d >= firstNonFixed || d <= -firstNonFixed) return null;

  // Find a sufficiently precise decimal representation of n.
  DecimalDecomposition decomposition = doubleToDecimal(
      d, DoubleConversionMode.fixed, requestedDigits: requestedDigits);

  return createDecimalRepresentation(decomposition, requestedDigits);
}

/// Converts the given double to an exponential representation.
///
/// If [requestedDigits] is -1 then the shortest representation in
/// exponential form is produced.
String convertToExponential(double d, int requestedDigits) {
  if (isSpecial(d)) return handleSpecialValues(d);

  if (requestedDigits < -1) return null;
  if (requestedDigits > maxExponentialDigits) return null;

  DecimalDecomposition decomposition;

  if (requestedDigits == -1) {
    decomposition = doubleToDecimal(d, DoubleConversionMode.shortest);
  } else {
    decomposition = doubleToDecimal(
        d, DoubleConversionMode.precision,
        requestedDigits: requestedDigits + 1);
    List<int> digits = decomposition.decimalDigits;
    int length = digits.length;
    // Fill the digits-list with 0s if there aren't enough digits.
    if (length < requestedDigits + 1) {
      digits.length = requestedDigits + 1;
      for (int i = length; i < digits.length; i++) {
        digits[i] = 0;
      }
    }
  }
  return createExponentialRepresentation(decomposition);
}

String convertToPrecision(double d, int precision) {
  if (isSpecial(d)) return handleSpecialValues(d);

  if (precision < minPrecisionDigits || precision > maxPrecisionDigits) {
    return null;
  }

  // Find a sufficiently precise decimal representation of n.
  DecimalDecomposition decomposition = doubleToDecimal(
      d, DoubleConversionMode.precision, requestedDigits: precision);

  int decimalPoint = decomposition.decimalPointPosition;
  List<int> digits = decomposition.decimalDigits;
  int length = digits.length;

  // The exponent if we print the number as x.xxeyyy. That is, with the
  // decimal point after the first digit.
  int exponent = decimalPoint - 1;

  if ((-decimalPoint + 1 > maxLeadingPaddingZeroesInPrecisionMode) ||
      (decimalPoint - precision > maxTrailingPaddingZeroesInPrecisionMode)) {
    // Fill buffer to contain 'precision' digits.
    // Usually the buffer is already at the correct length, but
    // [doubleToDecimal] is allowed to return fewer characters.
    if (length < precision) {
      digits.length = precision;
      for (int i = length; i < digits.length; i++) {
        digits[i] = 0;
      }
    }
    return createExponentialRepresentation(decomposition);
  }
  return  createDecimalRepresentation(
      decomposition, math.max(0, precision - decimalPoint));
}

String createDecimalRepresentation(
    DecimalDecomposition decomposition, int digitsAfterPoint) {

  // Adds the zero charcode [count] times to the list. If count is negative
  // nothing happens.
  void addZeroPadding(List<int> charCodes, int count) {
    for (int i = 0; i < count; i++) {
      charCodes.add(codes.zero);
    }
  }

  List<int> digits = decomposition.decimalDigits;
  int length = digits.length;
  int decimalPoint = decomposition.decimalPointPosition;
  bool isNegative = decomposition.isNegative;

  List<int> charCodes = <int>[];

  if (isNegative) charCodes.add(codes.hyphen);

  // Create a representation that is padded with zeros if needed.
  if (decimalPoint <= 0) {
    // "0.00000decimal_rep".
    charCodes.add(codes.zero);
    if (digitsAfterPoint > 0) {
      charCodes.add(codes.dot);
      addZeroPadding(charCodes, -decimalPoint);
      assert(length <= digitsAfterPoint - (-decimalPoint));
      for (int i = 0; i < digits.length; i++) {
        charCodes.add(codes.zero + digits[i]);
      }
      int remainingDigits = digitsAfterPoint - (-decimalPoint) - length;
      addZeroPadding(charCodes, remainingDigits);
    }
  } else if (decimalPoint >= length) {
    // "decimal_rep0000.00000" or "decimal_rep.0000"
    for (int i = 0; i < digits.length; i++) {
      charCodes.add(codes.zero + digits[i]);
    }
    addZeroPadding(charCodes, decimalPoint - length);
    if (digitsAfterPoint > 0) {
      charCodes.add(codes.dot);
      addZeroPadding(charCodes, digitsAfterPoint);
    }
  } else {
    // "decima.l_rep000"
    assert(digitsAfterPoint > 0);
    for (int i = 0; i < decimalPoint; i++) {
      charCodes.add(codes.zero + digits[i]);
    }
    charCodes.add(codes.dot);
    assert(length - decimalPoint <= digitsAfterPoint);
    for (int i = decimalPoint; i < digits.length; i++) {
      charCodes.add(codes.zero + digits[i]);
    }
    int remainingDigits = digitsAfterPoint - (length - decimalPoint);
    addZeroPadding(charCodes, remainingDigits);
  }
  return new String.fromCharCodes(charCodes);
}

// Constructs an exponential representation (i.e. 1.234e56).
// The given exponent assumes a decimal point after the first decimal digit.
String createExponentialRepresentation(DecimalDecomposition decomposition) {
  List<int> digits = decomposition.decimalDigits;
  int decimalPoint = decomposition.decimalPointPosition;
  bool isNegative = decomposition.isNegative;
  int exponent = decimalPoint - 1;
  int length = digits.length;
  assert(length != 0);

  List<int> charCodes = <int>[];

  if (isNegative) charCodes.add(codes.hyphen);

  for (int i = 0; i < digits.length; i++) {
    if (i == 1) charCodes.add(codes.dot);
    charCodes.add(codes.zero + digits[i]);
  }
  charCodes.add(codes.e);
  if (exponent < 0) {
    charCodes.add(codes.hyphen);
    exponent = -exponent;
  } else {
    charCodes.add(codes.plus);
  }
  if (exponent == 0) {
    charCodes.add(codes.zero);
  } else {
    assert(exponent < 1e4);
    // For simplicity add the digits in reverse order first and revert them
    // once we have all the digits.
    int firstExponentDigitPosition = charCodes.length;
    while (exponent > 0) {
      int exponentDigit = exponent.remainder(10);
      charCodes.add(codes.zero + exponentDigit);
      exponent = exponent ~/ 10;
    }
    int lastExponentDigitPosition = charCodes.length - 1;
    // Reverse the order of the exponent digits.
    while (firstExponentDigitPosition < lastExponentDigitPosition) {
      int tmp = charCodes[firstExponentDigitPosition];
      charCodes[firstExponentDigitPosition] =
          charCodes[lastExponentDigitPosition];
      charCodes[lastExponentDigitPosition] = tmp;
      firstExponentDigitPosition++;
      lastExponentDigitPosition--;
    }
  }
  return new String.fromCharCodes(charCodes);
}

// When calling convertToFixed with a double > 10^maxFixedDigitsBeforePoint
// or a requestedDigits parameter > maxFixedDigitsAfterPoint then the
// function returns null.
const int maxFixedDigitsBeforePoint = 60;
const int maxFixedDigitsAfterPoint = 60;

// If the decimal decomposition of a number has exponents in the range
// [lowExponentForDecimal] to [highExponentForDecimal], low inclusive, and high
// exclusive, then the number is printed as a decimal representation. Otherwise,
// it is printed as exponential.
const int lowExponentForDecimal = -6;
const int highExponentForDecimal = 21;

// When calling convertToExponential with a requestedDigits
// parameter > maxExponentialDigits then the function returns null.
const int maxExponentialDigits = 120;

// When calling convertToPrecision with a requestedDigits
// parameter < minPrecisionDigits or requestedDigits > maxPrecisionDigits
// then the function returns null.
const int minPrecisionDigits = 1;
const int maxPrecisionDigits = 120;

// [convertToPrecision] prefers to emit numbers in decimal mode. It allows
// [maxLeadingPaddingZeroesInPrecisionMode] and
// [maxTrailingPaddingZeroesInPrecisionMode] to make this happen.
// For example: convertToPrecision(0.1, 2) wants to precision digits. Since
// leading digits don't count as precision digits, it has to either return:
// 0.10 or 1.0e-1.
// Here, there would be one leading padding zero (and no trailing one), so
// [convertToPrecision] emits the decimal representation.
const int maxLeadingPaddingZeroesInPrecisionMode = 6;
const int maxTrailingPaddingZeroesInPrecisionMode = 0;

enum DoubleConversionMode {
  shortest,
  fixed,
  exponential,
  precision
}

bool isSpecial(double d) {
  // TODO(floitsch): we could use a fast bit-check for that.
  return d.isNaN || d.isInfinite;
}

String handleSpecialValues(double d) {
  if (d.isNaN) return "NaN";
  assert(d.isInfinite);
  return d < 0.0 ? "-Infinity" : "Infinity";
}

class DecimalDecomposition {
  final List<int> decimalDigits;
  final int decimalPointPosition;
  final bool isNegative;  // True if negative.

  DecimalDecomposition(
      this.decimalDigits, this.decimalPointPosition, this.isNegative);

  String toString() {
    return "$decimalDigits $decimalPointPosition $isNegative";
  }
}

DecimalDecomposition doubleToDecimal(
    double d, DoubleConversionMode mode, {int requestedDigits: -1}) {
  assert(!isSpecial(d));
  assert(mode == DoubleConversionMode.shortest || requestedDigits >= 0);

  bool isNegative = d.isNegative;  // Also handles -0.0.
  if (isNegative) d = -d;

  if (mode == DoubleConversionMode.precision && requestedDigits == 0) {
    return new DecimalDecomposition([], 0, isNegative);
  }

  if (d == 0) {
    return new DecimalDecomposition([0], 1, isNegative);
  }

  DoubleProperties properties = new DoubleProperties(d);
  int significand = properties.significand;
  int exponent = properties.exponent;
  bool isLowerBoundaryCloser = properties.isLowerBoundaryCloser;

  bool needBoundaryDeltas = mode == DoubleConversionMode.shortest;
  bool isEven = significand.isEven;
  int normalizedExponent = normalizeExponent(significand, exponent);
  int estimatedPower =
      estimatePower(normalizedExponent, DoubleProperties.significandSize);

  // Shortcut for fixed mode.
  // The requested digits correspond to the digits after the point. If the
  // number is way too small, then there is no need in trying to any digits.
  if (mode == DoubleConversionMode.fixed &&
      -estimatedPower - 1 > requestedDigits) {
    return new DecimalDecomposition([], -requestedDigits, isNegative);
  }

  ScaledStartValues startValues = computeScaledStartValues(
      significand, exponent, isLowerBoundaryCloser, estimatedPower,
      needBoundaryDeltas);

  int numerator = startValues.numerator;
  int denominator = startValues.denominator;
  int deltaMinus = startValues.deltaMinus;
  int deltaPlus = startValues.deltaPlus;

  // Multiply numerator/denominator so that its values lies in the
  // range 1-10. That is after a call to this function we have:
  //    1 <= (numerator + deltaPlus) /denominator < 10.
  // Let numerator the input before modification and numerator' the argument
  // after modification, then the output-parameter decimal_point is such that
  //  numerator / denominator * 10^estimatedPower ==
  //    numerator' / denominator' * 10^(decimalPoint - 1)
  // In some cases estimatedPower was too low, and this is already the case. We
  // then simply adjust the power so that 10^(k-1) <= v < 10^k (with k ==
  // estimatedPower) but do not touch the numerator or denominator.
  // Otherwise multiply the numerator and the deltas by 10.
  bool isInRange = isEven
      ? numerator + deltaPlus >= denominator
      : numerator + deltaPlus > denominator;

  int decimalPoint;
  if (isInRange) {
    decimalPoint = estimatedPower + 1;
  } else {
    decimalPoint = estimatedPower;
    numerator *= 10;
    deltaMinus *= 10;
    deltaPlus *= 10;
  }

  switch (mode) {
    case DoubleConversionMode.shortest:
      return doubleToDecimalShortest(
          numerator, denominator, deltaMinus, deltaPlus, decimalPoint,
          isNegative, isEven);
    case DoubleConversionMode.fixed:
      return doubleToDecimalFixed(
          requestedDigits, numerator, denominator, decimalPoint, isNegative);
    case DoubleConversionMode.precision:
      return doubleToDecimalPrecision(
          requestedDigits, numerator, denominator, decimalPoint, isNegative);
    default:
      throw "shouldn't reach $mode";
  }
}

// The procedure starts generating digits from the left to the right and stops
// when the generated digits yield the shortest decimal representation of v. A
// decimal representation of v is a number lying closer to v than to any other
// double, so it converts to v when read.
//
// This is true if d, the decimal representation, is between m- and m+, the
// upper and lower boundaries. d must be strictly between them if !is_even.
//           m- := (numerator - delta_minus) / denominator
//           m+ := (numerator + delta_plus) / denominator
//
// Precondition: 0 <= (numerator+delta_plus) / denominator < 10.
//   If 1 <= (numerator+delta_plus) / denominator < 10 then no leading 0 digit
//   will be produced. This should be the standard precondition.
DecimalDecomposition doubleToDecimalShortest(
    int numerator,
    int denominator,
    int deltaMinus,
    int deltaPlus,
    int decimalPoint,
    bool isNegative,
    bool isEven) {
  List<int> digits = <int>[];
  while(true) {
    // Usually the following two operations can be done in one go.
    int digit = numerator ~/ denominator;
    numerator = numerator.remainder(denominator);
    assert(0 <= digit && digit <= 9);
    digits.add(digit);

    // Can we stop already?
    // If the remainder of the division is less than the distance to the lower
    // boundary we can stop. In this case we simply round down (discarding the
    // remainder).
    // Similarly we test if we can round up (using the upper boundary).
    bool inDeltaRoomMinus = isEven
        ? numerator <= deltaMinus
        : numerator < deltaMinus;
    bool inDeltaRoomPlus = isEven
        ? numerator + deltaPlus >= denominator
        : numerator + deltaPlus > denominator;

    if (!inDeltaRoomMinus && !inDeltaRoomPlus) {
      // Prepare for next iteration.
      numerator *= 10;
      deltaMinus *= 10;
      deltaPlus *= 10;
    } else if (inDeltaRoomMinus && inDeltaRoomPlus) {
      // Let's see if denominator > 2*numerator
      // If yes, then the next digit would be < 5 and we can round down.
      int diff = (denominator - numerator - numerator);
      if (diff < 0 ||
          (diff == 0 && (digits.last & 1) == 1)) {
        // Either the next digit would be > 5, or it would be exactly 5.
        // In the halfway case we need to round to even. If the last digit is
        // not even, then we need to round up (for the half-way case).
        // Note that the last digit could not be a '9' as otherwise the whole
        // loop would have stopped earlier.
        assert(digits.last != 9);
        digits[digits.length - 1]++;
      } else {
        // The next digit would be < 5, or == 5, but we already round to even.
        // Round down (i.e. do nothing).
      }
      break;
    } else if (inDeltaRoomMinus) {
      // Round down. Nothing to do.
      break;
    } else if (inDeltaRoomPlus) {
      // Round up.
      // The last digit couldn't be '9' since this would have stopped the loop
      // earlier.
      assert(digits.last != 9);
      digits[digits.length - 1]++;
      break;
    }
  }
  return new DecimalDecomposition(digits, decimalPoint, isNegative);
}

// Generates 'requested_digits' after the decimal point. It might omit
// trailing '0's. If the input number is too small then no digits at all are
// generated (ex.: 2 fixed digits for 0.00001).
//
// Input verifies:  1 <= (numerator + delta) / denominator < 10.
DecimalDecomposition doubleToDecimalFixed(
    int requestedDigits,
    int numerator,
    int denominator,
    int decimalPoint,
    bool isNegative) {
  List<int> digits = <int>[];

  // Note that we have to look at more than just the requested_digits, since
  // a number could be rounded up. Example: v=0.5 with requested_digits=0.
  // Even though the power of v equals 0 we can't just stop here.
  if (-decimalPoint > requestedDigits) {
    // The number is definitively too small.
    // Ex: 0.001 with requested_digits == 1.
    // Set decimal-point to -requested_digits. This is what Gay does.
    // Note that it should not have any effect anyways since the string is
    // empty.
    return new DecimalDecomposition(<int>[], -requestedDigits, isNegative);
  }
  if (-decimalPoint == requestedDigits) {
    // We only need to verify if the number rounds down or up.
    // Ex: 0.04 and 0.06 with requested_digits == 1.
    // Initially the fraction lies in range (1, 10]. Multiply the denominator
    // by 10 so that we can compare more easily.
    if (numerator + numerator >= denominator * 10) {
      // If the fraction is >= 0.5 then we have to include the rounded
      // digit.
      return new DecimalDecomposition(<int>[1], decimalPoint + 1, isNegative);
    } else {
      // All zero.
      return new DecimalDecomposition(<int>[], -requestedDigits, isNegative);
    }
  }
  // The requested digits correspond to the digits after the point.
  // The variable 'needed_digits' includes the digits before the point.
  int neededDigits = decimalPoint + requestedDigits;
  return doubleToDecimalCounted(
      neededDigits, numerator, denominator, decimalPoint, isNegative);
}

DecimalDecomposition doubleToDecimalPrecision(
    int requestedDigits,
    int numerator,
    int denominator,
    int decimalPoint,
    bool isNegative) {
  return doubleToDecimalCounted(
      requestedDigits, numerator, denominator, decimalPoint, isNegative);
}

// Let v = numerator / denominator < 10.
// Then we generate 'count' digits of d = x.xxxxx... (without the decimal point)
// from left to right. Once 'count' digits have been produced we decide wether
// to round up or down. Remainders of exactly .5 round upwards. Numbers such
// as 9.999999 propagate a carry all the way, and change the
// exponent (decimalPoint), when rounding upwards.
DecimalDecomposition doubleToDecimalCounted(
    int count,
    int numerator,
    int denominator,
    int decimalPoint,
    bool isNegative) {
  assert(count >= 0);
  if (count == 0) {
    return new DecimalDecomposition(<int>[], decimalPoint, isNegative);
  }
  List<int> digits = <int>[];
  // Note that we don't compute the last digit.
  for (int i = 0; i < count - 1; ++i) {
    // Usually the following two operations can be done in one go.
    int digit = numerator ~/ denominator;
    numerator = numerator.remainder(denominator);
    assert(0 <= digit && digit <= 9);
    digits.add(digit);
    // Prepare for next iteration.
    numerator *= 10;
  }
  // Compute the last digit.
  // Adjust the last digit if necessary.
  int digit = numerator ~/ denominator;
  numerator = numerator.remainder(denominator);
  if (numerator + numerator >= denominator) {
    digit++;
  }
  digits.add(digit);
  // If the last digit is now a 10, propagate the carry back forward.
  for (int i = count - 1; i > 0; --i) {
    if (digits[i] == 10) {
      digits[i] = 0;
      digits[i - 1]++;
    }
  }
  if (digits[0] == 10) {
    digits[0] = 1;
    decimalPoint++;
  }
  return new DecimalDecomposition(digits, decimalPoint, isNegative);
}

class ScaledStartValues {
  // The following numbers have to be bignums.
  // The smallest double equals 4e-324. In that case the denominator needs fewer
  // than 324*4 binary digits. The maximum double is 1.7976931348623157e308
  // which needs fewer than 308*4 binary digits.
  // The bignum implementation must therefore be able to support ~1500 bits.
  final int numerator;
  final int denominator;
  final int deltaMinus;
  final int deltaPlus;

  ScaledStartValues(
      this.numerator, this.denominator, this.deltaMinus, this.deltaPlus);
}

// Let v = significand * 2^exponent.
// Computes v / 10^estimated_power exactly, as a ratio of two bignums, numerator
// and denominator.
//
// The initial start values consist of:
//  - a scaled numerator: s.t. numerator/denominator == v / 10^estimated_power.
//  - a scaled (common) denominator.
//  optionally:
//  - v - m-: the distance to the lower boundary.
//  - m+ - v: the distance to the upper boundary.
//
// v, m+, m-, and therefore v - m- and m+ - v all share the same denominator.
//
// Let ep == estimatedPower, then the returned values will satisfy:
//  v / 10^ep = numerator / denominator.
//  v's boundarys m- and m+:
//    m- / 10^ep == v / 10^ep - delta_minus / denominator
//    m+ / 10^ep == v / 10^ep + delta_plus / denominator
//  Or in other words:
//    m- == v - delta_minus * 10^ep / denominator;
//    m+ == v + delta_plus * 10^ep / denominator;
//
// Since 10^(k-1) <= v < 10^k    (with k == estimated_power)
//  or       10^k <= v < 10^(k+1)
//  we then have 0.1 <= numerator/denominator < 1
//           or    1 <= numerator/denominator < 10
ScaledStartValues computeScaledStartValues(
    int significand,
    int exponent,
    bool isLowerBoundaryCloser,
    int estimatedPower,
    bool needBoundaryDeltas) {
  ScaledStartValues result;
  if (exponent >= 0) {
    result = computeInitialScaledStartValuesPositiveExponent(
        significand, exponent, estimatedPower, needBoundaryDeltas);
  } else if (estimatedPower >= 0) {
    result = computeInitialScaledStartValuesNegativeExponentPositivePower(
        significand, exponent, estimatedPower, needBoundaryDeltas);
  } else {
    result = computeInitialScaledStartValuesNegativeExponentNegativePower(
        significand, exponent, estimatedPower, needBoundaryDeltas);
  }

  if (needBoundaryDeltas && isLowerBoundaryCloser) {
    // The lower boundary is closer at half the distance of "normal" numbers.
    // Increase the common denominator and adapt all but the delta_minus.
    return new ScaledStartValues(
      result.numerator * 2,
      result.denominator * 2,
      result.deltaMinus,
      result.deltaPlus * 2);
  }
  return result;
}

ScaledStartValues computeInitialScaledStartValuesPositiveExponent(
    int significand,
    int exponent,
    int estimatedPower,
    bool needBoundaryDeltas) {
  // A positive exponent implies a positive power.
  assert(estimatedPower >= 0);
  // Since the estimated_power is positive we simply multiply the denominator
  // by 10^estimatedPower.

  // numerator = v.
  int numerator = significand << exponent;
  // denominator = 10^estimatedPower.
  int denominator = math.pow(10, estimatedPower);

  int deltaPlus = 0;
  int deltaMinus = 0;
  if (needBoundaryDeltas) {
    // Introduce a common denominator so that the deltas to the boundaries are
    // integers.
    denominator <<= 1;
    numerator <<= 1;
    // Let v = f * 2^e, then m+ - v = 1/2 * 2^e; With the common
    // denominator (of 2) deltaPlus equals 2^e.
    deltaPlus = 1 << exponent;
    // Same for deltaMinus. The adjustments if f == 2^p-1 are done later.
    deltaMinus = 1 << exponent;
  }

  return new ScaledStartValues(numerator, denominator, deltaMinus, deltaPlus);
}

ScaledStartValues computeInitialScaledStartValuesNegativeExponentPositivePower(
    int significand,
    int exponent,
    int estimatedPower,
    bool needBoundaryDeltas) {
  // v = f * 2^e with e < 0, and with estimatedPower >= 0.
  // This means that e is close to 0 (have a look at how estimatedPower is
  // computed).

  // numerator = significand
  //  since v = significand * 2^exponent this is equivalent to
  //  numerator = v * / 2^-exponent
  int numerator= significand;
  // denominator = 10^estimatedPower * 2^-exponent (with exponent < 0)
  int tenToEstimatedPower = math.pow(10, estimatedPower);
  int denominator = tenToEstimatedPower << -exponent;

  int deltaPlus = 0;
  int deltaMinus = 0;
  if (needBoundaryDeltas) {
    // Introduce a common denominator so that the deltas to the boundaries are
    // integers.
    denominator <<= 1;
    numerator <<= 1;
    // Let v = f * 2^e, then m+ - v = 1/2 * 2^e; With the common
    // denominator (of 2) deltaPlus equals 2^e.
    // Given that the denominator already includes v's exponent the distance
    // to the boundaries is simply 1.
    deltaPlus = 1;
    // Same for deltaMinus. The adjustments if f == 2^p-1 are done later.
    deltaMinus = 1;
  }

  return new ScaledStartValues(numerator, denominator, deltaMinus, deltaPlus);
}

ScaledStartValues computeInitialScaledStartValuesNegativeExponentNegativePower(
    int significand,
    int exponent,
    int estimatedPower,
    bool needBoundaryDeltas) {
  // Instead of multiplying the denominator with 10^estimatedPower we
  // multiply all values (numerator and deltas) by 10^-estimatedPower.

  int powerTen = math.pow(10, -estimatedPower);

  int deltaPlus = 0;
  int deltaMinus = 0;

  if (needBoundaryDeltas) {
    // Since power_ten == numerator we must make a copy of 10^estimatedPower
    // before we complete the computation of the numerator.
    // delta_plus = delta_minus = 10^estimatedPower
    deltaPlus = powerTen;
    deltaMinus = powerTen;
  }

  // numerator = significand * 10^-estimatedPower
  //  since v = significand * 2^exponent this is equivalent to
  // numerator = v * 10^-estimatedPower * 2^-exponent.
  int numerator = powerTen * significand;

  // denominator = 2^-exponent with exponent < 0.
  int denominator = 1 << -exponent;

  if (needBoundaryDeltas) {
    // Introduce a common denominator so that the deltas to the boundaries are
    // integers.
    numerator <<= 1;
    denominator <<= 1;
    // With this shift the boundaries have their correct value, since
    // deltaPlus = 10^-estimatedPower, and
    // deltaMinus = 10^-estimatedPower.
    // These assignments have been done earlier.
    // The adjustments if f == 2^p-1 (lower boundary is closer) are done later.
  }

  return new ScaledStartValues(numerator, denominator, deltaMinus, deltaPlus);
}

/// Returns an estimation of k such that 10^(k-1) <= v < 10^k where
/// v = f * 2^exponent and 2^(p-1) <= f < 2^p.
/// v is hence a normalized double with the given exponent. The output is an
/// approximation for the exponent of the decimal approximation .digits * 10^k.
///
/// The result might undershoot by 1 in which case 10^k <= v < 10^k+1.
/// Note: this property holds for v's upper boundary m+ too.
///    10^k <= m+ < 10^k+1.
///   (see explanation below).
///
/// Examples:
///
///    estimatePower(0, 53)   => 16
///    estimatePower(-52, 53) => 0
///
/// Note: e >= 0 => estimatePower(e) > 0. No similar claim can be made for e<0.
int estimatePower(int exponent, int p) {
  // This function estimates log10 of v where v = f*2^e (with e == exponent).
  // Note that 10^floor(log10(v)) <= v, but v <= 10^ceil(log10(v)).
  // Note that f is bounded by its container size. For example, for a double
  // the significand has 53 bits. Thus p = 53 and 2^(p-1) <= f < 2^p.
  //
  // Given that log10(v) == log2(v)/log2(10) and e+(len(f)-1) is quite close
  // to log2(v) the function is simplified to (e+(len(f)-1)/log2(10)).
  // The computed number undershoots by less than 1/log2(10) =~= 0.301.
  //
  // Explanation for v's boundary m+: the computation takes advantage of
  // the fact that 2^(p-1) <= f < 2^p. Boundaries still satisfy this requirement
  // (even for denormals where the delta can be much more important).

  int undershootBias = 2000;
  // 1/ln(10) =~= 0.30102999566398_11952137388947244930267681898814621085
  // We use an approximated version: 0x134413509f79f * 2^-50
  //   =~= 0.30102999566398_036535019855364225804805755615234375
  //
  // We start by computing (e - p - 1) * 0x134413509f79f, and will deal with the
  // 2^-50 at the end (just before ceiling).
  //
  // Remember that we guarantee that the value is either correct or undershoots.
  // For positive numbers that's guaranteed, because the approximation is
  // smaller than the actual value. However, for
  // negative numbers this would lead to a potential overshoot. We correct for
  // this by decrementing the result of the multiplication by 2000.
  // We chose 2000 because abs(e - p - 1) < 2000 (for reasonable p). For
  // positive numbers we therefore undershoot even more, but it's still tiny
  // since we will multiply by 2^-50 in the next step.
  //
  // We combine the multiplication of 2^-50 with the ceiling operation.
  // Shifting out the least significant 50 bits is equivalent to flooring the
  // result. By decrementing the result before-hand (if only slightly) means
  // that adding 1 to the result of the shift is equivalent to ceiling the
  // full multiplication.
  assert((exponent + p - 1).abs() < undershootBias);
  return (((exponent + p - 1) * 0x134413509f79f - undershootBias) >> 50) + 1;
}

// Normalizes the exponent for denormalized values.
int normalizeExponent(int significand, int exponent) {
  assert(significand > 0);
  int currentBit = DoubleProperties.hiddenBit;
  while ((significand & currentBit) == 0) {
    currentBit >>= 1;
    exponent--;
  }
  return exponent;
}
