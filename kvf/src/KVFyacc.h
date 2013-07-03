/* A Bison parser, made by GNU Bison 2.3.  */

/* Skeleton interface for Bison's Yacc-like parsers in C

   Copyright (C) 1984, 1989, 1990, 2000, 2001, 2002, 2003, 2004, 2005, 2006
   Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.  */

/* As a special exception, you may create a larger work that contains
   part or all of the Bison parser skeleton and distribute that work
   under terms of your choice, so long as that work isn't itself a
   parser generator using the skeleton or a modified version thereof
   as a parser skeleton.  Alternatively, if you modify or redistribute
   the parser skeleton itself, you may (at your option) remove this
   special exception, which will cause the skeleton and the resulting
   Bison output files to be licensed under the GNU General Public
   License without this special exception.

   This special exception was added by the Free Software Foundation in
   version 2.2 of Bison.  */

/* Tokens.  */
#ifndef YYTOKENTYPE
# define YYTOKENTYPE
   /* Put the tokens into the symbol table, so that GDB and other debuggers
      know about them.  */
   enum yytokentype {
     T_CHARACTER_LITERAL = 258,
     T_INTEGER_LITERAL = 259,
     T_FLOATING_PT_LITERAL = 260,
     T_STRING_LITERAL = 261,
     T_IDENTIFIER = 262,
     T_LEFT_CURLY_BRACKET = 263,
     T_LEFT_PARANTHESIS = 264,
     T_LEFT_SQUARE_BRACKET = 265,
     T_RIGHT_CURLY_BRACKET = 266,
     T_RIGHT_PARANTHESIS = 267,
     T_RIGHT_SQUARE_BRACKET = 268,
     T_SEMI_COLON = 269,
     T_EQUALS = 270,
     T_COMMA = 271
   };
#endif
/* Tokens.  */
#define T_CHARACTER_LITERAL 258
#define T_INTEGER_LITERAL 259
#define T_FLOATING_PT_LITERAL 260
#define T_STRING_LITERAL 261
#define T_IDENTIFIER 262
#define T_LEFT_CURLY_BRACKET 263
#define T_LEFT_PARANTHESIS 264
#define T_LEFT_SQUARE_BRACKET 265
#define T_RIGHT_CURLY_BRACKET 266
#define T_RIGHT_PARANTHESIS 267
#define T_RIGHT_SQUARE_BRACKET 268
#define T_SEMI_COLON 269
#define T_EQUALS 270
#define T_COMMA 271




#if ! defined YYSTYPE && ! defined YYSTYPE_IS_DECLARED
typedef union YYSTYPE
#line 19 "KVFyacc.y"
{
  char         *ident;
  int            _int;
  char         *_str;
  char           _char;
  double         _float;
}
/* Line 1489 of yacc.c.  */
#line 89 "KVFyacc.h"
	YYSTYPE;
# define yystype YYSTYPE /* obsolescent; will be withdrawn */
# define YYSTYPE_IS_DECLARED 1
# define YYSTYPE_IS_TRIVIAL 1
#endif

extern YYSTYPE KVFyylval;

