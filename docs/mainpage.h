/**
 * @file   mainpage.h
 * @author Abtin Rahimian <arahimian@acm.org>
 * @date   Tue Jan 19 11:22:23 2010
 *
 * @brief  The main page for the Ves3D package manual.
 */


/*!@mainpage Introduction

  Some general info.

  This manual is divided in the following sections:
  - @subpage CodingConv
  - @subpage QuickStart
  - @subpage Developers
  - @subpage filelist
*/


////////////////////////////////////////////////////////////

/*!@page CodingConv Coding Conventions

  @section namingRules Naming Rules:
  We adopt a simple naming convention as following, for a more
  detailed explanation you can see the website <a
  href="http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml">
  here</a>.

  - <b> Type Names:</b> Type names start with a capital letter and
  have a capital letter for each new word, with no underscores:
  MyClass, MyEnum.  The names of all types—classes, structs, typedefs,
  and enums—have the same naming convention.<i> No underscores.</i>

  - <b>Variable Names:</b> Variable names are all lowercase, with
  underscores between words. Class member variables have trailing
  underscores. For instance, <code>my_local_variable,
  my_member_variable_.</code>

  - <b>Class Data Members:</b> Data members, also called instance
  variables or member variables, are treated like regular variable
  names, but always end with a <i>trailing underscore.</i>

  - <b>Struct Variables:</b> Data members in structs should be named
  like regular variables without the trailing underscores that data
  members in classes have.

  - <b>Global Variables:</b> They should be rare in any case, but if
  you use one, prefix it with <code>the_</code> marker to easily
  distinguish it from local variables.

  - <b>Constant Names:</b> All compile-time constants, whether they
  are declared locally, globally, or as part of a class, follow a
  slightly different naming convention from other variables. Use a
  <code>k_</code> followed by words, like <code>const int k_days_in_a_week =
  7</code>;

  - <b>Function Names:</b> Regular functions have mixed case;
    accessors and mutators match the name of the variable:
    MyExcitingFunction(), MyExcitingMethod(),
    my_exciting_member_variable(), set_my_exciting_member_variable().

    Regular functions should start with a capital letter and have a
    capital letter for each new word. No underscore, for example <code>
    AddTableEntry()</code> or <code> DeleteUrl()</code>.

    Accessors and mutators (get and set functions) should match the name
    of the variable they are getting and setting.

  - <b>Namespace Names:</b> Namespace names are all lower-case, and
  based on project names and possibly their directory structure.
*/

////////////////////////////////////////////////////////////

/*!@page QuickStart Quick Start
  This page introduces the user to the topic.
  Now you can proceed to the \ref advanced "advanced section".

  @section requiredPack Required Packages
  @section Install Installation
 */

////////////////////////////////////////////////////////////

/*!@page Developers Developers
  This page is for advanced users.
  Make sure you have first read \ref intro "the introduction".
 */
