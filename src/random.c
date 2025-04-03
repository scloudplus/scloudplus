#include <stdlib.h>

#if defined(_WIN32) || defined(_WIN64)

#include <windows.h>
#include <wincrypt.h>

#else
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#endif

int randombytes(unsigned char *buffer, unsigned int size) {
    if (buffer == NULL || size == 0) {
        return -1;  // 无效参数
    }

#if defined(_WIN32) || defined(_WIN64)
    /* Windows 实现 (兼容 MinGW/MSVC) */
    HCRYPTPROV hCryptProv = 0;

    // 获取加密服务提供程序句柄
    if (!CryptAcquireContext(&hCryptProv, NULL, NULL, PROV_RSA_FULL, CRYPT_VERIFYCONTEXT)) {
        return -2;  // 获取上下文失败
    }

    // 生成随机数据
    if (!CryptGenRandom(hCryptProv, size, buffer)) {
        CryptReleaseContext(hCryptProv, 0);
        return -3;  // 随机生成失败
    }

    // 释放上下文
    CryptReleaseContext(hCryptProv, 0);
    return 0;

#else
    /* UNIX-like 系统 (macOS/Linux) */
    int fd = open("/dev/urandom", O_RDONLY);
    if (fd == -1) {
        return -4;  // 打开设备失败
    }

    ssize_t bytes_read = 0;
    while (bytes_read < (ssize_t)size) {
        ssize_t result = read(fd, buffer + bytes_read, size - bytes_read);
        if (result < 0) {
            if (errno == EINTR) continue;  // 被信号中断则重试
            close(fd);
            return -5;  // 读取失败
        }
        bytes_read += result;
    }

    close(fd);
    return 0;
#endif
}