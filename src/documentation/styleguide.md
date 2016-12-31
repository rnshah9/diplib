# Style guide for contributors {#styleguide}

When contributing to *DIPlib*, please follow the style and layout of other files.
Specifically:

## Programming style:

- Everything should be declared within the `dip` namespace, or a sub-namespace. The
  exception is functionality that interfaces *DIPlib* with other libraries or software,
  which should be defined in their own namespaces (e.g. the `dml` namespace for the
  *DIPlib--MATLAB* interface).

- Don't put `dip::` in front of every identifyer within the library code, but do always
  do so for `dip::sint` and `dip::uint`, as they might be confused with types commonly
  present in the base namespace or as preprocessor macros.

- Do explicitly state the namespace for identifyers from every other library, including
  `std::`. This makes it easier to find references to specific types or functions.

- All functions local to a translation unit must be declared `static` to prevent
  name space pollution. This also prevents them for being exported out of the library.

- Prefer `using` over `typedef`.

- Use `struct` for classes without any private members.

- Option parameters to high-level functions (those that should be available in interfaces
  to other languages such as *MATLAB*) should be strings or string arrays, which are easier
  to translate to scripted languages.

- Option parameters to low-level functions (those that are meant to be called only from
  C++ code) should be defined through `DIP_DECLARE_OPTIONS`/`DIP_DEFINE_OPTION` or as
  `enum class`, preferably in `dip::Option::` or another sub-namespace. These are simpler
  and more efficient than strings.

- Don't use `bool` as a function parameter, prefer "yes"/"no" or "on"/"off" strings in
  high-level functions, and `enum class` with two options defined in `dip::Option::`
  namespace for low-level functions.

- Multiple return values are preferably combined in a `struct`, rather than a `std::tuple`
  or similar, as a `struct` has named members and is easier to use. Output should rarely
  be put into the function's argument list, with the exception of images
  (see \ref design_function_signatures).

## Naming conventions:

- Use camel case for variable, function and class names. Variable names start with
  a lowercase letter, function and class names start with an uppercase letter. Don't
  use underscores except for in a few special cases. Private class member variables
  end in an underscore. Enumerator constants are in all uppercase letters.

- Setter member functions start with `Set`. But getter member functions do not start
  with `Get`. Query functions that return a boolean start with `Is`. Member functions
  that do something have a name that resembles a verb: `Forge`, `Convert`, `PermuteDimensions`.

- The exception is in classes such as `dip::DimensionArray`, which is meant to emulate
  the `std::vector` class, and therefore follows the naming convention of the C++ Standard
  Library. Also, most classes define a `swap` operator that needs to be named as in
  the C++ Standard Library to be useful.

- Use all uppercase letters for preprocessor macros. Separate words with an underscore.
  Macros always start with `DIP_`, or `DIP__` if it is an internal macro not meant to
  be used outside the scope of the file in which it is defined. Include guards also
  start with `DIP_`, and end with `_H`; the part in the middle is an all-uppercase
  version of the file name.

- File names are in all lowercase and use underscores between words. There's no need
  to shorten names to 8 characters, so don't make the names cryptic.

## Formatting:

- All loops and conditional statements should be surrounded by braces, even if they
  are only one statement long.

- Indents are three spaces, don't use tab characters. Continuation indents are double
  regular indents, as are the indents for class definitions. Namespace scope is not
  indented.

- The keyword `const` comes after the type name it modifies: `dip::Image const& img`.

- Braces and brackets have spaces on the inside, not the outside.