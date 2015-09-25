library atod;

import 'dart:math' as math;
import 'char_codes.dart' as codes;
import 'float.dart';

// Every ascii character can be downcased by ORing it with 0x20.
const int downCaseBit = 0x20;

// Maximum number of significant digits in decimal representation.
// The longest possible double in decimal representation is
// (2^53 - 1) * 2 ^ -1074 that is (2 ^ 53 - 1) * 5 ^ 1074 / 10 ^ 1074
// (768 digits). If we parse a number whose first digits are equal to a
// mean of 2 adjacent doubles (that could have up to 769 digits) the result
// must be rounded to the bigger one unless the tail consists of zeros, so
// we don't need to preserve all the digits.
const int maxSignificantDigits = 772;

// Since the significand and the exponent can cancel each other out, we have to
// be prepared to accept decimal exponents that are out of normal ranges.
// In theory we could allow more, but the string would need to be so big that
// it becomes unreasonable.
const int maxDecimalExponent = 0x3FFFFFFF;

double atod(String str) {
  if (str == "NaN") return double.NAN;
  if (str == "Infinity") return double.INFINITY;
  if (str == "-Infinity") return double.NEGATIVE_INFINITY;
  if (str == "") return null;

  int startPos = 0;
  bool isNegative = false;
  if (str.codeUnitAt(0) == codes.hyphen) {
    isNegative = true;
    startPos = 1;
  }

  List<int> digits = <int>[];
  digits.length = str.length;

  int decimalPointPos = -1;
  int decimalSignificandLength = 0;
  int exponentPos = -1;
  int digitsPos = 0;
  for (int i = startPos; i < str.length; i++) {
    int codeUnit = str.codeUnitAt(i);
    if (codeUnit == codes.dot) {
      if (decimalPointPos != -1) return null;
      decimalPointPos = i;
      continue;
    } else if (codeUnit | downCaseBit == codes.e) {
      exponentPos = i;
      break;
    } else if (codes.zero <= codeUnit && codeUnit <= codes.nine) {
      digits[digitsPos++] = codeUnit - codes.zero;
    } else {
      return null;
    }
  }

  digits.length = digitsPos;

  bool exponentIsNegative = false;
  int exponent = 0;
  if (exponentPos != -1) {
    exponentPos++; // Jump over 'e'.

    // We need at least one digit.
    if (exponentPos == str.length) return null;
    int possibleSignUnit = str.codeUnitAt(exponentPos);
    if (possibleSignUnit == codes.hyphen) {
      exponentIsNegative = true;
      exponentPos++;
    } else if (possibleSignUnit == codes.plus) {
      // Ignore sign, but advance.
      exponentPos++;
    }
    // We still need at least one digit.
    if (exponentPos == str.length) return null;

    for (int i = exponentPos; i < str.length; i++) {
      int codeUnit = str.codeUnitAt(i);
      if (codes.zero <= codeUnit && codeUnit <= codes.nine) {
        exponent = exponent * 10 + (codeUnit - codes.zero);
      }
      if (exponent > maxDecimalExponent) return null;
    }
    if (exponentIsNegative) exponent = -exponent;
  }

  // Add the decimal-point position to the exponent. We won't need that
  // position afterwards.
  if (decimalPointPos != -1) {
    exponent -= (digits.length - decimalPointPos);
  }

  int nonZeroStart = 0;
  {
    int i;
    for (i = 0; i < digits.length; i++) {
      if (digits[i] != 0) break;
    }
    if (i == digits.length) return isNegative ? -0.0 : 0.0;
    nonZeroStart = i;
  }

  int nonZeroEnd = digits.length - 1;
  // We are guaranteed to encounter a non-zero digit. Therefore we don't need
  // to test against 0.
  while (digits[nonZeroEnd] == 0) nonZeroEnd--;

  // We can remove leading zeroes without needing to adjust the exponent, but
  // we need to pay attention when removing trailing zeroes.

  for (int i = nonZeroStart; i <= nonZeroEnd; i++) {
    digits[i - nonZeroStart] = digits[i];
  }
  exponent += digits.length - (nonZeroEnd + 1);
  digits.length -= nonZeroStart + digits.length - (nonZeroEnd + 1);

  // The smallest double is 4.9406564584124654*10^-324 and the
  // biggest double is 1.7976931348623157*10^308.
  // If the input is far from these numbers we can shortcut the computation.
  if (digits.length + exponent > 310) {
    return isNegative ? double.NEGATIVE_INFINITY : double.INFINITY;
  } else if (digits.length + exponent < -326) {
    return isNegative ? -0.0 : 0.0;
  }

  // We know we only need a limited number of digits from the decimal
  // representation. After some time we only need to know if a digit is
  // zero or not.
  if (digits.length > maxSignificantDigits) {
    // The trailing digits are non-zero (since we trimmed the array)
    int oldLength = digits.length;
    digits.length = maxSignificantDigits + 1;
    // In case the digits were 5,0,0,0,....,0,7 and we just cut of the trailing
    // non-zero digit. Add a non-zero digit as last kept digit.
    digits[maxSignificantDigits] = 1;
    // Adjust the exponent.
    exponent += oldLength - digits.length;
  }

  // Now to the complicated part of transforming the digits into a double.

  // The following is a shortcut. It can be removed (but it should trigger for
  // many numbers).
  {
    double result = convertDigitsFastPath(digits, exponent);
    if (result != null) return isNegative ? -result : result;
  }

  double result = convertDigits(digits, exponent);
  return isNegative ? -result : result;
}


// 2^53 = 9007199254740992.
// Any integer with at most 15 decimal digits will hence fit into a double
// (which has a 53bit significand) without loss of precision.
const int maxExactDoubleIntegerDecimalDigits = 15;

// If size is important this table can be replaced by a call to
// `math.pow(2.0, i)`.
const List<double> exactPowersOfTen = const <double>[
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

int convertDigitsToInt(List<int> digits) {
  int result = 0;
  for (int i = 0; i < digits.length; i++) {
    result = result * 10 + digits[i];
  }
  return result;
}

// If the digits and the 10^exponent fit into a double without loss of
// precision we can compute the double directly with floating-point operations.
double convertDigitsFastPath(List<int> digits, int exponent) {
  // If we have few digits, then they fit into a double without loss of
  // precision.
  if (digits.length >= maxExactDoubleIntegerDecimalDigits) return null;

  if (0 <= exponent && exponent < exactPowersOfTen.length) {
    // 10^exponent fits into a double.
    return convertDigitsToInt(digits) * exactPowersOfTen[exponent];
  }
  if (-(exactPowersOfTen.length) < exponent && exponent < 0) {
    // 10^-exponent fits into a double.
    return convertDigitsToInt(digits) / exactPowersOfTen[-exponent];
  }
  // The exponent is outside the exact-powers-of-ten range, but maybe the
  // significand is short enough that we can adjust it.
  int remainingDigits = maxExactDoubleIntegerDecimalDigits - digits.length;
  if (0 <= exponent && exponent - remainingDigits < exactPowersOfTen.length) {
    // The trimmed string was short and we can multiply it with
    // 10^remaining_digits. As a result the remaining exponent now fits
    // into a double too.
    int digitsAsInteger = convertDigitsToInt(digits);
    double shiftedDigits = digitsAsInteger * exactPowersOfTen[remainingDigits];
    int shiftedExponent = exponent - remainingDigits;
    return shiftedDigits * exactPowersOfTen[shiftedExponent];
  }
  // Shortcut didn't work.
  return null;
}

double convertDigits(List<int> digits, int exponent) {
  // Compute an approximation first, and then adjust it in a loop.

  // 17 digits fit into a 63 bit integer but is still more precise than what a
  // double can hold.
  int approximationDigitCount = digits.length < 17 ? digits.length : 17;
  int decimalSignificandApproximation = 0;
  for (int i = 0; i < approximationDigitCount; i++) {
    decimalSignificandApproximation =
        decimalSignificandApproximation * 10 + digits[i];
  }
  int adjustedExponent = exponent + (digits.length - approximationDigitCount);

  // The candidate works with imprecise approximations:
  //   the significant has been potentially truncated, and math.pow for the
  //   given values might not be the most precise approximation.
  //   Finally, the result of the multiplication or division can also introduce
  //   an error.
  // However, the error should be relatively small, only a few ulps (units in
  // the last place).
  double candidate;
  if (adjustedExponent >= 0) {
    candidate = decimalSignificandApproximation.toDouble() *
        math.pow(10.0, adjustedExponent);
  } else if (adjustedExponent <= -308) {
    // If the exponent is too small, math.pow(10, -adjustedExponent) would yield
    // infinity. Do the computation in two steps, to avoid this.
    candidate = decimalSignificandApproximation.toDouble() * 1e-308 /
        math.pow(10.0, -adjustedExponent - 308);
  } else {
    assert(-308 < adjustedExponent && adjustedExponent < 0);
    candidate = decimalSignificandApproximation.toDouble() /
        math.pow(10.0, -adjustedExponent);
  }
  // It's easier if the candidate is never 0.
  if (candidate == 0.0) candidate = DoubleProperties.minPositiveNonZeroDouble;

  int decimalSignificand = 0;
  for (int i = 0; i < digits.length; i++) {
    // TODO(floitsch): we could avoid repetitive bigint operations by computing
    // the digits in bigger chunks.
    decimalSignificand = decimalSignificand * 10 + digits[i];
  }

  int tenToPositiveDecimalExponent = math.pow(10, exponent.abs());

  // Iteratively move towards the correct double.
  //
  // Let v = digits * 10^exponent.
  // We have to check if v is within the range of the candidate. This comparison
  // happens with bigints, so we avoid expensive computations (like divisions).
  assert(candidate != 0.0);
  while (true) {
    // If the candidate is zero here, then it was decremented in a previous
    // iteration. Therefore it wasn't in the boundary of the next double.
    if (candidate == 0.0) return candidate;

    // Note that we decompose infinity as if it had a value. For IEEE
    // floating-point numbers, infinity is represented by the biggest
    // representation. For the purpose of the computations below we use the
    // "normal" value of that representation, instead of treating it like
    // infinity.
    DoubleProperties properties = new DoubleProperties(candidate);
    int candidateSignificand = properties.significand;
    int candidateExponent = properties.exponent;

    // To avoid divisions we scale the inputs (v and the candidate) with a
    // positive scale-factor:
    // scaledV = v * scaleFactor and
    // scaledCandidate = candidate * scaleFactor.
    // For example, for a negative exponent and a positive candidateExponent,
    // the scaleFactor would be 10^-exponent.
    int scaledV, scaledCandidate;
    if (exponent >= 0 && candidateExponent >= 0) {
      scaledV = decimalSignificand * tenToPositiveDecimalExponent;
      scaledCandidate = candidateSignificand * (1 << candidateExponent);
    } else if (exponent >= 0 && candidateExponent < 0) {
      scaledV = decimalSignificand * tenToPositiveDecimalExponent *
          (1 << -candidateExponent);
      scaledCandidate = candidateSignificand;
    } else if (exponent < 0 && candidateExponent >= 0) {
      scaledV = decimalSignificand;
      scaledCandidate = candidateSignificand * tenToPositiveDecimalExponent *
          (1 << candidateExponent);
    } else {
      assert(exponent < 0 && candidateExponent < 0);
      scaledV = decimalSignificand * (1 << -candidateExponent);
      scaledCandidate = candidateSignificand * tenToPositiveDecimalExponent;
    }

    int diff = scaledV - scaledCandidate;
    if (diff == 0) return candidate;
    if (diff > 0 && candidate == double.INFINITY) return candidate;

    if (diff > 0) {
      // v is bigger than the candidate, but maybe it's within the boundaries
      // of candidate. For a given normalized floating-point number d=m*2^k, the
      // next-bigger number d2 = (m+1)*2^k. The half-way point is the boundary
      // between the two: boundary = (d2-d)/2. Resolving the equation we find
      // boundary = d + 2^(k-1).
      // We want to check if v-candidate < 2^(candidateExponent-1). We want to
      // avoid computing 2^(candidateExponent-1), so we start by multiplying
      // both sides by 2:
      // (v-candidate)*2 < 2^candidateExponent. We don't have
      // 2^candidateExponent either, but it's one of the factors of the
      // candidate: candidate = candidateSignificand * 2^candidateExponent. By
      // multiplying with candidateSignificand we can replace the right-hand
      // side with candidate:
      //   (v-candidate)*2*candidateSignificand < candidate.
      // Note that we can only do this if the candidateSignificand is not 0.
      assert(candidateSignificand != 0);
      int scaledBoundaryDiff = diff * candidateSignificand * 2;
      if (scaledBoundaryDiff < scaledCandidate) {
        // The candidate is within the boundary.
        return candidate;
      } else if (scaledBoundaryDiff == scaledCandidate) {
        // The candidate lies exactly on the boundary. We round to even.
        if (candidateSignificand.isEven) return candidate;
        return properties.nextDouble;
      } else {
        // Try the next double.
        candidate = properties.nextDouble;
      }
    } else {
      // v is smaller than the candidate.
      // Similar to above we have to check if v lies within the boundaries of
      // the candidate. However, this time we have to pay more attention,
      // because the lower boundary might be closer, if the candidate's
      // significand is 1<<p (where p is the double's significand size).
      // Example (in decimal, with significand size equal to 5):
      //   candidate = 10000*10^3. Then the next lower floating point number
      //   would be 99999*10^2, at distance 1*2^2 (instead of 1*2^3).
      int scaledBoundaryDiff;
      // For doubles the hidden bit is equal to 1<<p.
      if (candidateSignificand == DoubleProperties.hiddenBit) {
        // Boundary is at distance 2^(candidateExponent-2). Otherwise everything
        // works as above.
        scaledBoundaryDiff = -diff * candidateSignificand * 4;
      } else {
        scaledBoundaryDiff = -diff * candidateSignificand * 2;
      }
      if (scaledBoundaryDiff < scaledCandidate) {
        // The candidate is within the boundaries.
        return candidate;
      } else if (scaledBoundaryDiff == scaledCandidate) {
        // The candidate lies exactly on the boundary. We round to even.
        if (candidateSignificand.isEven) return candidate;
        return properties.previousDouble;
      } else {
        // Try the next double.
        candidate = properties.previousDouble;
      }
    }
  }
}
