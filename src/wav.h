/*!
 * @brief wav.cのヘッダファイル
 *
 * @author  koturn 0;
 * @date    2013 06/01
 * @file    wav.h
 * @version 0.5
 */
#ifndef WAV_H
#define WAV_H

#include "compatibility.h"


typedef struct {
  char            riff_chunk_id[4];     // "RIFF"
  unsigned int    riff_chunk_size;
  char            file_format_type[4];  // "WAVE"

  char            fmt_chunk_id[4];      // "fmt "
  unsigned int    fmt_chunk_size;       // 16
  unsigned short  wave_format_type;     // 1 : PCM
  unsigned short  channel;              // 1 : Mono, 2 : Stereo
  unsigned int    samples_per_sec;      // Sampling frequency
  unsigned int    bytes_per_sec;        // block_size * sample_per_sec
  unsigned short  block_size;           // bit_per_sample * channel / 8
  unsigned short  bits_per_sample;      // digitalize

  char            data_chunk_id[4];     // "data"
  unsigned int    data_chunk_size;      // sample_length * data_chunk_size;
  double         *s_data1;
  double         *s_data2;

  unsigned int    length;               // number of samples
} WAVE_DATA;


WAVE_DATA *read_wave_file(const char *restrict file_name);
void       free_wave_data(WAVE_DATA *restrict wd);
void       print_wave_header(const WAVE_DATA *restrict wd);

WAVE_DATA *make_skelton_wave_data(void);
int set_wave_data(
    WAVE_DATA      *restrict wd,
    const double   *restrict s_data1,
    const double   *restrict s_data2,
    unsigned int             N,
    unsigned int             samples_per_sec,
    unsigned short           channel,
    unsigned short           bits_per_sample);
void write_wave_file(const WAVE_DATA *restrict wd, const char *restrict file_name);


static const double RATE_8BIT  = 128.0;
static const double RATE_16BIT = 32768.0;

//! ステレオ、モノラルを示す列挙体
enum SOUND_TYPE {
  MONORAUL = 1,  //!< モノラル
  STEREO   = 2   //!< ステレオ
};

//! 量子化ビット数を示す列挙体
enum BITS_PER_SEC {
  BIT_8  = 8,  //!< 量子化ビット数は8bit
  BIT_16 = 16  //!< 量子化ビット数は16bit
};


#endif
