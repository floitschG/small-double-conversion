import 'dtoa.dart';

main(List<String> args) {
  print(convertToShortest(1.5));
  print(convertToShortest(0.1));
  print(convertToShortest(1e20));
  print(convertToShortest(1e19));
  print(convertToShortest(1e21));
  print(convertToShortest(3.141));
  print(convertToShortest(1.0));

  print(convertToFixed(1.5, 2));
  print(convertToFixed(0.1, 2));
  print(convertToFixed(1e20, 2));
  print(convertToFixed(1e19, 2));
  print(convertToFixed(1e21, 2));
  print(convertToFixed(3.141, 2));
  print(convertToFixed(1.0, 2));

  print("precision");
  print(convertToPrecision(1.5, 2));
  print(convertToPrecision(0.1, 2));
  print(convertToPrecision(1e20, 2));
  print(convertToPrecision(1e19, 2));
  print(convertToPrecision(1e21, 2));
  print(convertToPrecision(3.141, 2));
  print(convertToPrecision(1.0, 2));
}
