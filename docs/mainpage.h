/**
 * @file   mainpage.h
 * @author Abtin Rahimian <abtin@romario>
 * @date   Tue Jan 19 11:22:23 2010
 * 
 * @brief  The main page for the Ves3D package manual.
 */


/*!@mainpage A simple manual
  
  Some general info.
  
  This manual is divided in the following sections:
  - @subpage CodingConv
  - @subpage Client
  - @subpage Developer
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
  accessors and mutators match the name of the variable. Functions
  should start with a capital letter and have a capital letter for
  each new word. No underscore, for example <code> AddTableEntry()</code>
  or <code> DeleteUrl()</code>.

  - <b>Namespace Names:</b> Namespace names are all lower-case, and
  based on project names and possibly their directory structure.
  
  @section emacs Personalizing Emacs
  
  I use emacs as the text editor and the following is the part in my
  <code>.emacs</code> file that corresponds to presonalization for C++:


  <code>
  (setq c-hungry-delete-key t)<br /> 
  (add-hook 'c-mode-common-hook 'flyspell-prog-mode)<br /> 

  ;; C++ mode style, space only as tab<br /> 
  (require 'cc-mode)<br /> 
  (defun my-build-tab-stop-list (width)<br /> 
  (let ((num-tab-stops (/ 80 width))<br /> 
  (counter 1)<br /> 
  (ls nil))<br /> 
  (while (<= counter num-tab-stops)<br /> 
  (setq ls (cons (* width counter) ls))<br /> 
  (setq counter (1+ counter)))<br /> 
  (set (make-local-variable 'tab-stop-list) (nreverse ls))))<br /> 
  
  (defun my-c-mode-common-hook ()<br /> 
  (setq tab-width 4)<br /> 
  (my-build-tab-stop-list tab-width)<br /> 
  (setq c-basic-offset tab-width)<br /> 
  (setq indent-tabs-mode nil)) ;; force only spaces for indentation<br /> 
  (add-hook 'c-mode-common-hook 'my-c-mode-common-hook)<br /> 

  ;; get the name of the symbol with C-c C-o<br /> 
  (c-set-offset 'substatement-open 0)<br /> 
  (c-set-offset 'case-label '+)<br /> 
  (c-set-offset 'access-label '-2)<br /> 
  (c-set-offset 'arglist-cont-nonempty '+)<br /> 
  (c-set-offset 'arglist-intro '+)<br /> 
  (c-set-offset 'topmost-intro-cont 0)<br /> 
  (c-set-offset 'topmost-intro 0)<br /> 

  ;; Doxygen<br /> 
  (add-hook 'c-mode-common-hook 'doxymacs-mode) <br /> 
  (defun my-doxymacs-font-lock-hook ()<br /> 
  (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))<br /> 
      (doxymacs-font-lock)))<br /> 
  (add-hook 'font-lock-mode-hook 'my-doxymacs-font-lock-hook)<br /> 

  </code>
*/

////////////////////////////////////////////////////////////

/*!@page Client Clients' Manual
  This page introduces the user to the topic.
  Now you can proceed to the \ref advanced "advanced section".
  
  @section requiredPack Required Packages
  @section Install Installation
 */

////////////////////////////////////////////////////////////

/*!@page Developer Developers' Documentation
  This page is for advanced users.
  Make sure you have first read \ref intro "the introduction".

  @section shTran Spherical Harmonics Transform
  @section diff Differentiation
 */
