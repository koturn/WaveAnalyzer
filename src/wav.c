/*!
 * @brief waveファイルを取り扱う関数を定義している
 *
 * @author  koturn 0;
 * @date    2013 06/01
 * @file    wav.c
 * @version 0.6
 */
#ifdef _MSC_VER
#  define _USE_MATH_DEFINES
#endif
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "compatibility.h"
#include "wav.h"

#define CALLOC(type, n)  ((type *)calloc((n), sizeof(type)))
#define FREE(ptr_p)      (free(*(ptr_p)), *(ptr_p) = NULL)


static void read_8bit_monoraul(WAVE_DATA *restrict wd, FILE *restrict fp);
static void read_16bit_monoraul(WAVE_DATA *restrict wd, FILE *restrict fp);
static void read_8bit_stereo(WAVE_DATA *restrict wd, FILE *restrict fp);
static void read_16bit_stereo(WAVE_DATA *restrict wd, FILE *restrict fp);
static void write_8bit_monoraul(const WAVE_DATA *restrict wd, FILE *restrict fp);
static void write_16bit_monoraul(const WAVE_DATA *restrict wd, FILE *restrict fp);
static void write_8bit_stereo(const WAVE_DATA *restrict wd, FILE *restrict fp);
static void write_16bit_stereo(const WAVE_DATA *restrict wd, FILE *restrict fp);

inline static void write_wave_header(const WAVE_DATA *restrict wd, FILE *restrict fp);
inline static void read_wave_header(WAVE_DATA *restrict WAVE_DATA, FILE *restrict fp);
inline static int  wave_err_check(const WAVE_DATA *restrict wd);
__attribute__((const)) inline static double clipping(double s, double max, double min);


static const double CLIPPING_MIN_8BIT  = 0.0;
static const double CLIPPING_MAX_8BIT  = 255.0;
static const double CLIPPING_MIN_16BIT = 0.0;
static const double CLIPPING_MAX_16BIT = 65535.0;


/*!
 * @brief waveファイルを読み込む
 *
 * サンプリング周波数やチャンネル数は自動判定する
 * @param [in] file_name  読み込むwaveファイル名
 * @return waveデータの構造体へのポインタ
 */
WAVE_DATA *read_wave_file(const char *restrict file_name) {
  FILE *restrict fp = fopen(file_name, "rb");
  WAVE_DATA *wd = CALLOC(WAVE_DATA, 1);
  read_wave_header(wd, fp);
  if (wd->channel == MONORAUL) {
    if (wd->bits_per_sample == 8) {
      read_8bit_monoraul(wd, fp);
    } else {
      read_16bit_monoraul(wd, fp);
    }
  } else {
    if (wd->bits_per_sample == 8) {
      read_8bit_stereo(wd, fp);
    } else {
      read_16bit_stereo(wd, fp);
    }
  }
  fclose(fp);
  return wd;
}


/*!
 * @brief waveファイルデータの構造体の領域を解放する
 * @param [in] wd  破棄する領域
 */
void free_wave_data(WAVE_DATA *restrict wd) {
  FREE(&wd->s_data1);
  FREE(&wd->s_data2);
  FREE(&wd);
}


/*!
 * @brief waveファイルのヘッダ情報を表示する
 * @param [in] wd  waveファイルデータ
 */
void print_wave_header(const WAVE_DATA *restrict wd) {
  printf("riff_chunk_id    (4byte):%12.4s (%s)\n", wd->riff_chunk_id, "constant");
  printf("riff_chunk_size  (4byte):%12u (%s)\n",   wd->riff_chunk_size, "36 + data_chunk_size");
  printf("file_format_type (4byte):%12.4s (%s)\n", wd->file_format_type, "constant");

  printf("fmt_chunk_id     (4byte):%12.4s (%s)\n", wd->fmt_chunk_id, "constant");
  printf("fmt_chunk_size   (4byte):%12u (%s)\n",   wd->fmt_chunk_size, "constant");
  printf("wave_format_type (2byte):%12u (%s)\n",   wd->wave_format_type, wd->wave_format_type == 1 ? "PCM" : "UNKNOWN");
  printf("channel          (2byte):%12u (%s)\n",   wd->channel, wd->channel == 1 ? "MONAURAL" : "STEREO");
  printf("samples_per_sec  (4byte):%12u\n",        wd->samples_per_sec);
  printf("bytes_per_sec    (4byte):%12u (%s)\n",   wd->bytes_per_sec, "block_size * samples_per_sec");
  printf("block_size       (2byte):%12u (%s)\n",   wd->block_size, "bits_per_sample * channel / 8");
  printf("bits_per_sample  (2byte):%12u\n",        wd->bits_per_sample);

  printf("data_chunk_id    (4byte):%12.4s (%s)\n", wd->data_chunk_id, "constant");
  printf("data_chunk_size  (4byte):%12u\n",        wd->data_chunk_size);
}


/*!
 * @brief 8bitモノラルのwaveデータを読み込む
 * @param [out] wd  格納先waveデータ構造体
 * @param [in]  fp  ファイルポインタ
 */
static void read_8bit_monoraul(WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned char data;  // 1byte
  unsigned int  len = wd->data_chunk_size;
  unsigned int  n;
  wd->s_data1 = CALLOC(double, len);
  wd->s_data2 = NULL;

  wd->length = len;
  for (n = 0; n < len; n++) {
    fread(&data, sizeof(data), 1, fp);                        // 音データの読み取り
    wd->s_data1[n] = ((double)data - RATE_8BIT) / RATE_8BIT;  // 音データを-1以上1未満の範囲に正規化する
  }
}


/*!
 * @brief 16bitモノラルのwaveデータを読み込む
 * @param [out] wd  格納先waveデータ構造体
 * @param [in]  fp  ファイルポインタ
 */
static void read_16bit_monoraul(WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int len = wd->data_chunk_size / 2;
  unsigned int n;
  short        data;  // 2byte
  wd->s_data1 = CALLOC(double, len);
  wd->s_data2 = NULL;

  wd->length = len;
  for (n = 0; n < len; n++) {
    fread(&data, sizeof(data), 1, fp);           // 音データの読み取り
    wd->s_data1[n] = (double)data / RATE_16BIT;  // 音データを-1以上1未満の範囲に正規化する
  }
}


/*!
 * @brief 8bitステレオのwaveデータを読み込む
 * @param [out] wd  格納先waveデータ構造体
 * @param [in]  fp  ファイルポインタ
 */
static void read_8bit_stereo(WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned char data;  // 1byte
  unsigned int  len = wd->data_chunk_size / 2;
  unsigned int  n;
  wd->s_data1 = CALLOC(double, len);
  wd->s_data2 = CALLOC(double, len);

  wd->length = len;
  for (n = 0; n < len; n++) {
    fread(&data, sizeof(data), 1, fp);                        // 音データの読み取り
    wd->s_data1[n] = ((double)data - RATE_8BIT) / RATE_8BIT;  // 音データを-1以上1未満の範囲に正規化する
    fread(&data, sizeof(data), 1, fp);                        // 音データの読み取り
    wd->s_data2[n] = ((double)data - RATE_8BIT) / RATE_8BIT;  // 音データを-1以上1未満の範囲に正規化する
  }
}


/*!
 * @brief 16bitステレオのwaveデータを読み込む
 * @param [out] wd  格納先waveデータ構造体
 * @param [in]  fp  ファイルポインタ
 */
static void read_16bit_stereo(WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int len = wd->data_chunk_size / 4;
  unsigned int n;
  short        data;  // 2byte
  wd->s_data1 = CALLOC(double, len);
  wd->s_data2 = CALLOC(double, len);

  wd->length = len;
  for (n = 0; n < len; n++) {
    fread(&data, sizeof(data), 1, fp);           // 音データの読み取り
    wd->s_data1[n] = (double)data / RATE_16BIT;  // 音データを-1以上1未満の範囲に正規化する
    fread(&data, sizeof(data), 1, fp);           // 音データの読み取り
    wd->s_data2[n] = (double)data / RATE_16BIT;  // 音データを-1以上1未満の範囲に正規化する
  }
}


/*!
 * @brief waveファイルのヘッダ情報を読み取る
 * @param [out] wd  格納先waveデータ構造体
 * @param [in]  fp  ファイルポインタ
 */
inline static void read_wave_header(WAVE_DATA *restrict wd, FILE *restrict fp) {
  fread( wd->riff_chunk_id,    sizeof(char),           4, fp);
  fread(&wd->riff_chunk_size,  sizeof(unsigned int),   1, fp);
  fread( wd->file_format_type, sizeof(char),           4, fp);

  fread( wd->fmt_chunk_id,     sizeof(char),           4, fp);
  fread(&wd->fmt_chunk_size,   sizeof(unsigned int),   1, fp);
  fread(&wd->wave_format_type, sizeof(unsigned short), 1, fp);
  fread(&wd->channel,          sizeof(unsigned short), 1, fp);
  fread(&wd->samples_per_sec,  sizeof(unsigned int),   1, fp);
  fread(&wd->bytes_per_sec,    sizeof(unsigned int),   1, fp);
  fread(&wd->block_size,       sizeof(unsigned short), 1, fp);
  fread(&wd->bits_per_sample,  sizeof(unsigned short), 1, fp);

  fread( wd->data_chunk_id,    sizeof(char),           4, fp);
  fread(&wd->data_chunk_size,  sizeof(unsigned int),   1, fp);
  /**
  int err;
  if ((err = wave_err_check(wd)) != 0) {
    fprintf(stderr, "invalid wave file! : error = %d\n", err);
  }
  **/
}


/*!
 * @brief WAVWデータが正しいフォーマットかどうかチェックする
 * @param [in] wd  検査するWAVEデータの構造体へのポインタ
 * @return 異常がないなら0を、異常があればそれに応じた負値を返す
 */
inline static int wave_err_check(const WAVE_DATA *restrict wd) {
  if (strncmp(wd->riff_chunk_id, "RIFF", sizeof(char) * 4)) {
    return -1;
  } else if (wd->riff_chunk_size != wd->data_chunk_size + 36) {
    return -2;
  } else if (strncmp(wd->file_format_type, "WAVE", sizeof(char) * 4)) {
    return -3;
  } else if (strncmp(wd->fmt_chunk_id, "fmt ", sizeof(char) * 4)) {
    return -4;
  } else if (wd->fmt_chunk_size != 16) {
    return -5;
  } else if (wd->wave_format_type != 1) {
    return -6;
  } else if (wd->channel != 1 && wd->channel != 2) {
    return -7;
  } else if (wd->bytes_per_sec != wd->block_size * wd->samples_per_sec) {
    return -8;
  } else if (wd->block_size != wd->bits_per_sample * wd->channel / 8) {
    return -9;
  } else {
    return 0;
  }
}


/*!
 * @brief 空のWAVWデータ構造体を作成する
 * @return 作成したWAVEデータの構造体へのポインタ
 */
WAVE_DATA *make_skelton_wave_data(void) {
  WAVE_DATA *wd = CALLOC(WAVE_DATA, 1);
  strncpy(wd->riff_chunk_id,    "RIFF", sizeof(char) * 4);  // Constatn
  strncpy(wd->file_format_type, "WAVE", sizeof(char) * 4);  // Constatn

  strncpy(wd->fmt_chunk_id,     "fmt ", sizeof(char) * 4);  // "fmt "
  wd->fmt_chunk_size   = 16;                                // Constant
  wd->wave_format_type = 1;                                 // Constatnt  1 : PCM

  strncpy(wd->data_chunk_id,    "data", sizeof(char) * 4);  // Constant
  return wd;
}


/*!
 * @brief WAVWデータ構造体の各メンバに値をセットする
 * @param [out] wd               WAVEデータ構造体へのポインタ
 * @param [in]  s_data1          WAVEデータにセットする音声波形配列(Lチャンネル)
 * @param [in]  s_data2          WAVEデータにセットする音声波形配列(Rチャンネル/モノラルの場合、関係ない)
 * @param [in]  N                音声波形配列の要素数
 * @param [in]  samples_per_sec  音声波形配列のサンプリング周波数
 * @param [in]  channel          チャンネル数(1 : モノラル, 2 : ステレオ)
 * @param [in]  bits_per_sample  量子化ビット数
 * @return 正常終了時:0, 異常終了時:-1
 */
int set_wave_data(
    WAVE_DATA      *restrict wd,
    const double   *restrict s_data1,
    const double   *restrict s_data2,
    unsigned int             N,
    unsigned int             samples_per_sec,
    unsigned short           channel,
    unsigned short           bits_per_sample)
{
  wd->samples_per_sec = samples_per_sec;
  wd->channel         = channel;
  wd->bits_per_sample = bits_per_sample;
  wd->data_chunk_size = N * channel * (bits_per_sample / 8);

  wd->riff_chunk_size = 36 + wd->data_chunk_size;
  wd->block_size      = wd->bits_per_sample * wd->channel / 8;
  wd->bytes_per_sec   = wd->block_size * wd->samples_per_sec;

  if (s_data1 != NULL) {
    if ((wd->s_data1 = CALLOC(double, N)) == NULL) {
      return -1;
    }
    memcpy(wd->s_data1, s_data1, N);
  }
  if (s_data2 != NULL) {
    if ((wd->s_data2 = CALLOC(double, N)) == NULL) {
      return -1;
    }
    memcpy(wd->s_data2, s_data2, N);
  }
  return 0;
}


/*!
 * @brief WAVEデータをファイルに書き出す
 * @param [in] wd         ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] file_name  書き出すファイル名
 */
void write_wave_file(const WAVE_DATA *restrict wd, const char *restrict file_name) {
  FILE *restrict fp = fopen(file_name, "wb");

  write_wave_header(wd, fp);
  if (wd->channel == MONORAUL) {
    if (wd->bits_per_sample == 8) {
      write_8bit_monoraul(wd, fp);
    } else {
      write_16bit_monoraul(wd, fp);
    }
  } else {
    if (wd->bits_per_sample == 8) {
      write_8bit_stereo(wd, fp);
    } else {
      write_16bit_stereo(wd, fp);
    }
  }
  fclose(fp);
}


/*!
 * @brief WAVEファイルのヘッダ情報をファイルに出力する
 * @param [in] wd  ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] fp  書き出すファイルへのポインタ
 */
inline static void write_wave_header(const WAVE_DATA *restrict wd, FILE *restrict fp) {
  fwrite( wd->riff_chunk_id,    sizeof(char),           4, fp);
  fwrite(&wd->riff_chunk_size,  sizeof(unsigned int),   1, fp);
  fwrite( wd->file_format_type, sizeof(char),           4, fp);

  fwrite( wd->fmt_chunk_id,     sizeof(char),           4, fp);
  fwrite(&wd->fmt_chunk_size,   sizeof(unsigned int),   1, fp);
  fwrite(&wd->wave_format_type, sizeof(unsigned short), 1, fp);
  fwrite(&wd->channel,          sizeof(unsigned short), 1, fp);
  fwrite(&wd->samples_per_sec,  sizeof(unsigned int),   1, fp);
  fwrite(&wd->bytes_per_sec,    sizeof(unsigned int),   1, fp);
  fwrite(&wd->block_size,       sizeof(unsigned short), 1, fp);
  fwrite(&wd->bits_per_sample,  sizeof(unsigned short), 1, fp);

  fwrite( wd->data_chunk_id,    sizeof(char),           4, fp);
  fwrite(&wd->data_chunk_size,  sizeof(unsigned int),   1, fp);
}


/*!
 * @brief 8bitモノラルのWAVEデータをファイルに出力する
 * @param [in] wd  ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] fp  書き出すファイルへのポインタ
 */
static void write_8bit_monoraul(const WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int i;
  unsigned int len = wd->data_chunk_size;
  for (i = 0; i < len; i++) {
    unsigned char data;
    double s = (wd->s_data1[i] + 1.0) / 2.0 * (CLIPPING_MAX_8BIT + 1.0);
    s = clipping(s, CLIPPING_MAX_8BIT, CLIPPING_MIN_8BIT);
    data = (unsigned char)((int)(s + 0.5));  // 四捨五入
    fwrite(&data, sizeof(data), 1, fp);      // 音データの書き出し
  }
  if ((len % 2) == 1) {                  // 音データの長さが奇数のとき
    unsigned char data = 0;
    fwrite(&data, sizeof(data), 1, fp);  // 0パディング
  }
}


/*!
 * @brief 16bitモノラルのWAVEデータをファイルに出力する
 * @param [in] wd  ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] fp  書き出すファイルへのポインタ
 */
static void write_16bit_monoraul(const WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int i;
  unsigned int len = wd->data_chunk_size / 2;
  for (i = 0; i < len; i++) {
    short data;
    double s = (wd->s_data1[i] + 1.0) / 2.0 * (CLIPPING_MAX_16BIT + 1.0);
    s = clipping(s, CLIPPING_MAX_16BIT, CLIPPING_MIN_16BIT);
    data = (short)((int)(s + 0.5) - ((CLIPPING_MAX_16BIT + 1.0) / 2));  // 四捨五入とオフセットの調節
    fwrite(&data, sizeof(data), 1, fp);  // 音データの書き出し
  }
}


/*!
 * @brief 8bitステレオのWAVEデータをファイルに出力する
 * @param [in] wd  ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] fp  書き出すファイルへのポインタ
 */
static void write_8bit_stereo(const WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int i;
  unsigned int len = wd->data_chunk_size / 2;

  for (i = 0; i < len; i++) {
    unsigned char data;
    double sL, sR;
    sL = (wd->s_data1[i] + 1.0) / 2.0 * (CLIPPING_MAX_8BIT + 1.0);
    sL = clipping(sL, CLIPPING_MAX_8BIT, CLIPPING_MIN_8BIT);
    data = (unsigned char)((int)(sL + 0.5));  // 四捨五入 */
    fwrite(&data, sizeof(data), 1, fp);       // 音データ（Lチャンネル）の書き出し

    sR = (wd->s_data2[i] + 1.0) / 2.0 * (CLIPPING_MAX_8BIT + 1.0);
    sR = clipping(sR, CLIPPING_MAX_8BIT, CLIPPING_MIN_8BIT);
    data = (unsigned char)((int)(sR + 0.5));  // 四捨五入
    fwrite(&data, sizeof(data), 1, fp);       // 音データ（Rチャンネル）の書き出し
  }
}


/*!
 * @brief 16bitステレオのWAVEデータをファイルに出力する
 * @param [in] wd  ファイルに出力するWAVEデータ構造体へのポインタ
 * @param [in] fp  書き出すファイルへのポインタ
 */
static void write_16bit_stereo(const WAVE_DATA *restrict wd, FILE *restrict fp) {
  unsigned int i;
  unsigned int len = wd->data_chunk_size / 4;

  for (i = 0; i < len; i++) {
    short  data;
    double sL, sR;
    sL = (wd->s_data1[i] + 1.0) / 2.0 * (CLIPPING_MAX_16BIT + 1.0);
    sL = clipping(sL, CLIPPING_MAX_16BIT, CLIPPING_MIN_16BIT);
    data = (short)((int)(sL + 0.5) - ((CLIPPING_MAX_16BIT + 1.0) / 2));  // 四捨五入とオフセットの調節
    fwrite(&data, sizeof(data), 1, fp);  // 音データ（Lチャンネル）の書き出し

    sR = (wd->s_data2[i] + 1.0) / 2.0 * (CLIPPING_MAX_16BIT + 1.0);
    sR = clipping(sR, CLIPPING_MAX_16BIT, CLIPPING_MIN_16BIT);
    data = (short)((int)(sR + 0.5) - ((CLIPPING_MAX_16BIT + 1.0) / 2));  // 四捨五入とオフセットの調節
    fwrite(&data, sizeof(data), 1, fp);  // 音データ（Rチャンネル）の書き出し
  }
}


/*!
 * @brief クリッピングを行い、値を調節する
 * @param [in] s    クリッピング対象の値
 * @param [in] max  クリッピングの上限値
 * @param [in] min  クリッピングの下限値
 * @return クリッピング後の値
 */
inline static double clipping(double s, double max, double min) {
  return s > max ? max
                 : s < min ? min
                           : s;
}
