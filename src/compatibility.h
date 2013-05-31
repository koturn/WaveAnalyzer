/*!
 * @brief 異なる処理系間での互換性を保つマクロを提供する
 * @author  koturn 0;
 * @date    2013 06/01
 * @file    compatibility.h
 * @version 0.1
 */
#ifndef COMPATIBILITY_H
#define COMPATIBILITY_H


#ifdef inline
#  warning 'inline' was already defined.
#endif
#ifdef __inline
#  warning '__inline' was already defined.
#endif
#ifdef __inline__
#  warning '__inline__' was already defined.
#endif
#ifdef restrict
#  warning 'restrict' was already defined.
#endif
#ifdef __restrict
#  warning '__restrict' was already defined.
#endif
#ifdef __restrict__
#  warning '__restrict__' was already defined.
#endif
#ifdef __attribute__
#  warning 'attribute' was already defined.
#endif


// inline指定に関する修正
#ifndef __cplusplus
#  if defined(_MSC_VER)
#    define inline      __inline  // Visual C++では、inlineではなく、__inline
#    define __inline__  __inline  // __inline__も同様
#  elif !defined(__GNUC__) && __STDC_VERSION__ < 199901L
#    define inline                // inline指定が使えない処理系では、inline指定を消す
#    define __inline
#  endif
#endif

// restrict指定に関する修正
#if _MSC_VER >= 1400
#  define restrict      __restrict  // Visual C++のCコンパイラでは、restrictと__restrict__は
#  define __restrict__  __restrict  // 使えないが、__restrictは使える
#elif __cplusplus
#  define restrict      __restrict  // C++では、restrictは使えず、__restrict__は使える
#elif !defined(__STDC_VERSION__) || __STDC_VERSION__ < 199901L
#  if defined(__GNUC__)
#    define restrict    __restrict  // C99以前のgccではrestrictではなく、__restrict
#  else
#    define restrict                // restrictが使えない処理系では、restrict指定を消す
#    define __restrict
#    define __restrict__
#  endif
#endif

#ifndef __GNUC__
#  define __attribute__(attr)
#endif


#endif
