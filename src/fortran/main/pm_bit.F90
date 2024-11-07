!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!                                                                                                                            !!!!
!!!!    ParaMonte: Parallel Monte Carlo and Machine Learning Library.                                                           !!!!
!!!!                                                                                                                            !!!!
!!!!    Copyright (C) 2012-present, The Computational Data Science Lab                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!    This file is part of the ParaMonte library.                                                                             !!!!
!!!!                                                                                                                            !!!!
!!!!    LICENSE                                                                                                                 !!!!
!!!!                                                                                                                            !!!!
!!!!       https://github.com/cdslaborg/paramonte/blob/main/LICENSE.md                                                          !!!!
!!!!                                                                                                                            !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!>  \brief
!>  This module contains constants and procedures that are relevant to bit manipulation.<br>
!>
!>  \final
!>
!>  \author
!>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module pm_bit

    use pm_kind, only: SK, IK, LK

    !>  \brief
    !>  The scalar `integer` constant of intrinsic default kind that is `.true.` <b>if and only if</b>
    !>  the hardware upon which the ParaMonte library is built is a **big-endian** system.<br>
    !>
    !>  \details
    !>  In computing, **endianness** is the order or sequence of bytes of a word of digital data in computer memory.<br>
    !>  Endianness is primarily expressed as **big-endian (BE)** or **little-endian (LE)**:<br>
    !>  <ol>
    !>      <li>    A big-endian system stores the most significant byte of a word at
    !>              the smallest memory address and the least significant byte at the largest.<br>
    !>      <li>    A little-endian system, in contrast, stores the least-significant byte at the smallest address.<br>
    !>      <li>    **Bi-endianness** is a feature supported by numerous computer architectures that feature
    !>              switchable endianness in data fetches and stores or for instruction fetches.<br>
    !>      <li>    Other orderings are generically called **middle-endian** or **mixed-endian**.<br>
    !>  </ol>
    !>  The following graph schematically illustrates the difference between the big-endian and little-endian systems
    !>  (credit: [Wikipedia](https://en.wikipedia.org/wiki/Endianness)).<br>
    !>  <br>
    !>  \htmlonly
    !>      <img src="pm_bit@IS_BIG_ENDIAN.png" style="width:40%;">
    !>  \endhtmlonly
    !>  <br>
    !>
    !>  <b>Registers and endianness</b><br>
    !>  Endianness only makes sense when a multi-byte quantity is broken up into bytes to be stored at consecutive memory locations.<br>
    !>  However, if a 32-bit register stores a 32-bit value, it makes no sense to talk about endianness.<br>
    !>  The register is neither big-endian nor little-endian; it is just a register holding a 32-bit value.<br>
    !>  The **rightmost bit is the least significant bit**, and the **leftmost bit is the most significant bit**.<br>
    !>  Because it stores its most significant byte at the lowest memory address, the register is sometimes classified as big-endian.<br>
    !>
    !>  <b>Importance of endianness</b><br>
    !>  Endianness is the attribute of a system that indicates whether integers are represented from left to right or right to left.<br>
    !>  Endianness must be chosen every time a hardware or software architecture is designed.<br>
    !>  There is not much in the way of natural law to help decide, so implementations vary.<br>
    !>  All processors must be designated as either big-endian or little-endian.<br>
    !>  For example, the `80x86` (aka, `amd64`) processors from Intel® and their clones are little-endian.<br>
    !>  By contrast, the Sun SPARC, Motorola 68K, and the PowerPC® families are all big-endian.<br>
    !>  Why is endianness so important? Suppose you are storing integer values to a file,
    !>  and you send the file to a machine that uses the opposite endianness as it reads in the value.<br>
    !>  This causes problems because of endianness; The values will be read in reverse, that do not make sense.<br>
    !>  Endianness is also a big issue when sending numbers over the network.<br>
    !>  Again, if you send a value from a machine of one endianness to a machine of the opposite endianness, there will be problems.<br>
    !>  This is even worse over the network because you might not be able to determine the endianness of the machine that sent you the data.<br>
    !>
    !>  <b>Determining the endianness at run time</b><br>
    !>  One way to determine the endianness is to test the memory layout of a predefined constant.<br>
    !>  For example, the layout of a 32-bit integer variable with a value of `1` is `00 00 00 01` for big-endian and `01 00 00 00` for little-endian.<br>
    !>  By looking at the first byte of the constant, one can tell the endianness of the running platform and then take the appropriate action.<br>
    !>
    !>  <b>Advantages of different endianness systems</b><br>
    !>  There seems to be no significant advantage in using one method of endianness over the other.<br>
    !>  Both are still common and different architectures use them.<br>
    !>  Little-endian based processors (and their clones) are used in most personal computers and laptops.<br>
    !>  Therefore, the vast majority of desktop computers today are little-endian.<br>
    !>  Endian issues do not affect sequences that have single bytes, because **byte** is considered an atomic unit from a storage point of view.<br>
    !>  On the other hand, sequences based on multi-byte are affected by endianness and you need to take care while coding.<br>
    !>
    !>  \see
    !>  [Fortran Discourse](https://fortran-lang.discourse.group/t/moving-bits-question/4799/14)<br>
    !>  [How to write endian-independent code in C](https://developer.ibm.com/articles/au-endianc/)<br>
    !>
    !>  \final{IS_BIG_ENDIAN}
    !>
    !>  \author
    !>  \AmirShahmoradi, September 1, 2017, 11:35 PM, Institute for Computational Engineering and Sciences (ICES), The University of Texas Austin<br>
    logical(LK) , parameter :: IS_BIG_ENDIAN = iachar(c = transfer(source = 1, mold = "a")) == 0

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end module pm_bit ! LCOV_EXCL_LINE
