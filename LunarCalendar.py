import math
from math import pi, atan, sqrt, sin, cos


class LunarCalendar:
    zwz = True
    ctg = ['甲', '乙', '丙', '丁', '戊', '己', '庚', '辛', '壬', '癸']
    cdz = ['子', '丑', '寅', '卯', '辰', '巳', '午', '未', '申', '酉', '戌', '亥']
    jq = ['春分', '清明', '谷雨', '立夏', '小满', '芒种', '夏至', '小暑', '大暑', '立秋', '处暑', '白露', '秋分', '寒露', '霜降', '立冬', '小雪', '大雪',
          '冬至', '小寒', '大寒', '立春', '雨水', '惊蛰']
    synmonth = 29.530588853
    ptsa = [485, 203, 199, 182, 156, 136, 77, 74, 70, 58, 52, 50, 45, 44, 29, 18, 17, 16, 14, 12, 12, 12, 9, 8]
    ptsb = [324.96, 337.23, 342.08, 27.85, 73.14, 171.52, 222.54, 296.72, 243.58, 119.81, 297.17, 21.02, 247.54, 325.15,
            60.93, 155.12, 288.79, 198.04, 199.76, 95.39, 287.11, 320.81, 227.73, 15.45]
    ptsc = [1934.136, 32964.467, 20.186, 445267.112, 45036.886, 22518.443, 65928.934, 3034.906, 9037.513, 33718.147,
            150.678, 2281.226, 29929.562, 31555.956, 4443.417, 67555.328, 4562.452, 62894.029, 31436.921, 14577.848,
            31931.756, 34777.259, 1222.114, 16859.074]

    def VE(self, yy: int) -> float:
        if yy < -8000:
            return False
        if yy > 8001:
            return False
        if 1000 <= yy <= 8001:
            m = (yy - 2000) / 1000
            return 2451623.80984 + 365242.37404 * m + 0.05169 * m * m - 0.00411 * m * m * m - 0.00057 * m * m * m * m
        if -8000 <= yy < 1000:
            m = yy / 1000
            return 1721139.29189 + 365242.1374 * m + 0.06134 * m * m + 0.00111 * m * m * m - 0.00071 * m * m * m * m
        return False

    def Perturbation(self, jd) -> float:
        t = (jd - 2451545) / 36525
        s = 0
        for k in range(24):
            s = s + self.ptsa[k] * math.cos(self.ptsb[k] * 2 * math.pi / 360 + self.ptsc[k] * 2 * math.pi / 360 * t)
        w = 35999.373 * t - 2.47
        l = 1 + 0.0334 * math.cos(w * 2 * math.pi / 360) + 0.0007 * math.cos(2 * w * 2 * math.pi / 360)
        return 0.00001 * s / l

    def DeltaT(self, yy: int, mm: int) -> float:
        y = yy + (mm - 0.5) / 12
        if y <= -500:
            u = (y - 1820) / 100
            dt = (-20 + 32 * u ** 2)
        elif y < 500:
            u = y / 100
            dt = (
                    10583.6 - 1014.41 * u + 33.78311 * u ** 2 - 5.952053 * u ** 3 - 0.1798452 * u ** 4 + 0.022174192 * u ** 5 + 0.0090316521 * u ** 6)
        elif y < 1600:
            u = (y - 1000) / 100
            dt = (
                    1574.2 - 556.01 * u + 71.23472 * u ** 2 + 0.319781 * u ** 3 - 0.8503463 * u ** 4 - 0.005050998 * u ** 5 + 0.0083572073 * u ** 6)
        elif y < 1700:
            t = y - 1600
            dt = (120 - 0.9808 * t - 0.01532 * t ** 2 + t ** 3 / 7129)
        elif y < 1800:
            t = y - 1700
            dt = (8.83 + 0.1603 * t - 0.0059285 * t ** 2 + 0.00013336 * t ** 3 - t ** 4 / 1174000)
        elif y < 1860:
            t = y - 1800
            dt = (
                    13.72 - 0.332447 * t + 0.0068612 * t ** 2 + 0.0041116 * t ** 3 - 0.00037436 * t ** 4 + 0.0000121272 * t ** 5 - 0.0000001699 * t ** 6 + 0.000000000875 * t ** 7)
        elif y < 1900:
            t = y - 1860
            dt = (7.62 + 0.5737 * t - 0.251754 * t ** 2 + 0.01680668 * t ** 3 - 0.0004473624 * t ** 4 + t ** 5 / 233174)
        elif y < 1920:
            t = y - 1900
            dt = (-2.79 + 1.494119 * t - 0.0598939 * t ** 2 + 0.0061966 * t ** 3 - 0.000197 * t ** 4)
        elif y < 1941:
            t = y - 1920
            dt = (21.2 + 0.84493 * t - 0.0761 * t ** 2 + 0.0020936 * t ** 3)
        elif y < 1961:
            t = y - 1950
            dt = (29.07 + 0.407 * t - t ** 2 / 233 + t ** 3 / 2547)
        elif y < 1986:
            t = y - 1975
            dt = (45.45 + 1.067 * t - t * t / 260 - t * t * t / 718)
        elif y < 2005:
            t = y - 2000
            dt = (
                    63.86 + 0.3345 * t - 0.060374 * t * t + 0.0017275 * t * t * t + 0.000651814 * t * t * t * t + 0.00002373599 * t * t * t * t * t)
        elif y < 2050:
            t = y - 2000
            dt = (62.92 + 0.32217 * t + 0.005589 * t * t)
        elif y < 2150:
            u = (y - 1820) / 100
            dt = (-20 + 32 * u * u - 0.5628 * (2150 - y))
        else:
            u = (y - 1820) / 100
            dt = (-20 + 32 * u * u)
        if y < 1955 or y >= 2005:
            dt = dt - (0.000012932 * (y - 1955) * (y - 1955))
        return dt / 60

    def MeanJQJD(self, yy: int) -> list:
        jd = self.VE(yy)
        if not jd:
            return []

        ty = self.VE(yy + 1) - jd
        num = 24 + 2
        ath = 2 * math.pi / 24
        tx = (jd - 2451545) / 365250
        e = 0.0167086342 - 0.0004203654 * tx - 0.0000126734 * tx * tx + 0.0000001444 * tx * tx * tx - 0.0000000002 * tx * tx * tx * tx + 0.0000000003 * tx * tx * tx * tx * tx
        tt = yy / 1000
        vp = 111.25586939 - 17.0119934518333 * tt - 0.044091890166673 * tt * tt - 4.37356166661345E-04 * tt * tt * tt + 8.16716666602386E-06 * tt * tt * tt * tt
        rvp = vp * 2 * math.pi / 360
        peri = []
        for i in range(num):
            flag = 0
            th = ath * i + rvp
            if math.pi < th <= 3 * math.pi:
                th = 2 * math.pi - th
                flag = 1
            if th > 3 * math.pi:
                th = 4 * math.pi - th
                flag = 2
            f1 = 2 * math.atan((math.sqrt((1 - e) / (1 + e)) * math.tan(th / 2)))
            f2 = (e * math.sqrt(1 - e * e) * math.sin(th)) / (1 + e * math.cos(th))
            f = (f1 - f2) * ty / 2 / math.pi
            if flag == 1:
                f = ty - f
            if flag == 2:
                f = 2 * ty - f
            peri.append(f)
        jqjd = []
        for i in range(num):
            jqjd.append(jd + peri[i] - peri[0])

        return jqjd

    def GetAdjustedJQ(self, yy: int, start: int, end: int) -> list:
        if start < 0 or start > 25:
            return []
        if end < 0 or end > 25:
            return []

        jq = []

        jqjd = self.MeanJQJD(yy)
        for k, jd in enumerate(jqjd):
            if k < start:
                continue
            if k > end:
                continue
            ptb = self.Perturbation(jd)
            dt = self.DeltaT(yy, int((k + 1) / 2) + 3)
            jq.append(jd + ptb - dt / 60 / 24 + 1 / 3)

        return jq

    def GetPureJQsinceSpring(self, yy: int) -> list:
        jdpjq = []

        dj = self.GetAdjustedJQ(yy - 1, 19, 23)
        for k, v in enumerate(dj):
            if k < 0:
                continue
            if k > 4:
                continue
            if k % 2 == 0:
                continue
            jdpjq.append(dj[k])

        dj = self.GetAdjustedJQ(yy, 0, 25)
        for k, v in enumerate(dj):
            if k % 2 == 0:
                continue
            jdpjq.append(dj[k])

        return jdpjq

    def GetZQsinceWinterSolstice(self, yy: int) -> list:
        jdzq = []

        dj = self.GetAdjustedJQ(yy - 1, 18, 23)
        jdzq.append(dj[0])
        jdzq.append(dj[2])
        jdzq.append(dj[4])

        dj = self.GetAdjustedJQ(yy, 0, 23)
        for k, v in enumerate(dj):
            if k % 2 != 0:
                continue
            jdzq.append(dj[k])

        return jdzq

    def TrueNewMoon(self, k):
        jdt = 2451550.09765 + k * self.synmonth
        t = (jdt - 2451545) / 36525
        t2 = t ** 2
        t3 = t2 * t
        t4 = t3 * t
        pt = jdt + 0.0001337 * t2 - 0.00000015 * t3 + 0.00000000073 * t4
        m = 2.5534 + 29.10535669 * k - 0.0000218 * t2 - 0.00000011 * t3
        mprime = 201.5643 + 385.81693528 * k + 0.0107438 * t2 + 0.00001239 * t3 - 0.000000058 * t4
        f = 160.7108 + 390.67050274 * k - 0.0016341 * t2 - 0.00000227 * t3 + 0.000000011 * t4
        omega = 124.7746 - 1.5637558 * k + 0.0020691 * t2 + 0.00000215 * t3
        es = 1 - 0.002516 * t - 0.0000074 * t2

        apt1 = -0.4072 * math.sin(math.radians(mprime))
        apt1 += 0.17241 * es * math.sin(math.radians(m))
        apt1 += 0.01608 * math.sin(math.radians(2 * mprime))
        apt1 += 0.01039 * math.sin(math.radians(2 * f))
        apt1 += 0.00739 * es * math.sin(math.radians(mprime - m))
        apt1 -= 0.00514 * es * math.sin(math.radians(mprime + m))
        apt1 += 0.00208 * es * es * math.sin(math.radians(2 * m))
        apt1 -= 0.00111 * math.sin(math.radians(mprime - 2 * f))
        apt1 -= 0.00057 * math.sin(math.radians(mprime + 2 * f))
        apt1 += 0.00056 * es * math.sin(math.radians(2 * mprime + m))
        apt1 -= 0.00042 * math.sin(math.radians(3 * mprime))
        apt1 += 0.00042 * es * math.sin(math.radians(m + 2 * f))
        apt1 += 0.00038 * es * math.sin(math.radians(m - 2 * f))
        apt1 -= 0.00024 * es * math.sin(math.radians(2 * mprime - m))
        apt1 -= 0.00017 * math.sin(math.radians(omega))
        apt1 -= 0.00007 * math.sin(math.radians(mprime + 2 * m))
        apt1 += 0.00004 * math.sin(math.radians(2 * mprime - 2 * f))
        apt1 += 0.00004 * math.sin(math.radians(3 * m))
        apt1 += 0.00003 * math.sin(math.radians(mprime + m - 2 * f))
        apt1 += 0.00003 * math.sin(math.radians(2 * mprime + 2 * f))
        apt1 -= 0.00003 * math.sin(math.radians(mprime + m + 2 * f))
        apt1 += 0.00003 * math.sin(math.radians(mprime - m + 2 * f))
        apt1 -= 0.00002 * math.sin(math.radians(mprime - m - 2 * f))
        apt1 -= 0.00002 * math.sin(math.radians(3 * mprime + m))
        apt1 += 0.00002 * math.sin(math.radians(4 * mprime))

        apt2 = 0.000325 * math.sin(math.radians(299.77 + 0.107408 * k - 0.009173 * t2))
        apt2 += 0.000165 * math.sin(math.radians(251.88 + 0.016321 * k))
        apt2 += 0.000164 * math.sin(math.radians(251.83 + 26.651886 * k))
        apt2 += 0.000126 * math.sin(math.radians(349.42 + 36.412478 * k))
        apt2 += 0.00011 * math.sin(math.radians(84.66 + 18.206239 * k))
        apt2 += 0.000062 * math.sin(math.radians(141.74 + 53.303771 * k))
        apt2 += 0.00006 * math.sin(math.radians(207.14 + 2.453732 * k))
        apt2 += 0.000056 * math.sin(math.radians(154.84 + 7.30686 * k))
        apt2 += 0.000047 * math.sin(math.radians(34.52 + 27.261239 * k))
        apt2 += 0.000042 * math.sin(math.radians(207.19 + 0.121824 * k))
        apt2 += 0.00004 * math.sin(math.radians(291.34 + 1.844379 * k))
        apt2 += 0.000037 * math.sin(math.radians(161.72 + 24.198154 * k))
        apt2 += 0.000035 * math.sin(math.radians(239.56 + 25.513099 * k))
        apt2 += 0.000023 * math.sin(math.radians(331.55 + 3.592518 * k))
        return pt + apt1 + apt2

    def MeanNewMoon(self, jd):
        kn = math.floor((jd - 2451550.09765) / self.synmonth)
        jdt = 2451550.09765 + kn * self.synmonth
        t = (jdt - 2451545) / 36525
        thejd = jdt + 0.0001337 * t * t - 0.00000015 * t * t * t + 0.00000000073 * t * t * t * t
        return [kn, thejd]

    def Julian2Solar(self, jd):
        jd = float(jd)

        if jd >= 2299160.5:  # 1582年10月15日,此日起是儒略日历,之前是儒略历
            y4h = 146097
            init = 1721119.5
        else:
            y4h = 146100
            init = 1721117.5

        jdr = int(jd - init)
        yh = y4h / 4
        cen = int((jdr + 0.75) / yh)
        d = int(jdr + 0.75 - cen * yh)
        ywl = 1461 / 4
        jy = int((d + 0.75) / ywl)
        d = int(d + 0.75 - ywl * jy + 1)
        ml = 153 / 5
        mp = int((d - 0.5) / ml)
        d = int(d - 0.5 - 30.6 * mp + 1)
        y = int(100 * cen + jy)
        m = (mp + 2) % 12 + 1
        if m < 3:
            y = y + 1
        sd = int((jd + 0.5 - int(jd + 0.5)) * 24 * 60 * 60 + 0.00005)
        mt = int(sd / 60)
        ss = sd % 60
        hh = int(mt / 60)
        mt = mt % 60
        yy = int(y)
        mm = int(m)
        dd = int(d)

        return [yy, mm, dd, hh, mt, ss]

    def GetZQandSMandLunarMonthCode(self, yy: int):
        mc = []

        jdzq = self.GetZQsinceWinterSolstice(yy)
        jdnm = self.GetSMsinceWinterSolstice(yy, jdzq[0])
        yz = 0
        if int(jdzq[12] + 0.5) >= int(jdnm[13] + 0.5):
            for i in range(1, 15):
                if int(jdnm[i] + 0.5) > int(jdzq[i - 1 - yz] + 0.5) and int(jdnm[i + 1] + 0.5) <= int(
                        jdzq[i - yz] + 0.5):
                    mc.append(i - 0.5)
                    yz = 1
                else:
                    mc.append(i - yz)
        else:
            for i in range(13):
                mc.append(i)
            for i in range(13, 15):
                if int(jdnm[i] + 0.5) > int(jdzq[i - 1 - yz] + 0.5) and int(jdnm[i + 1] + 0.5) <= int(
                        jdzq[i - yz] + 0.5):
                    mc.append(i - 0.5)
                    yz = 1
                else:
                    mc.append(i - yz)
        return [jdzq, jdnm, mc]

    def GetSMsinceWinterSolstice(self, yy: int, jdws):
        tjd = []
        jd = self.Solar2Julian(yy - 1, 11, 1, 0, 0, 0)
        kn, thejd = self.MeanNewMoon(jd)
        for i in range(20):
            k = kn + i
            mjd = thejd + self.synmonth * i
            tjd.append(self.TrueNewMoon(k) + 1 / 3)
            tjd[i] = tjd[i] - self.DeltaT(yy, i - 1) / 1440

        j = 0
        while j <= 18 and math.floor(tjd[j] + 0.5) <= math.floor(jdws + 0.5):
            j += 1

        jdnm = []
        for k in range(16):
            jdnm.append(tjd[j - 1 + k])
        return jdnm

    def Solar2Julian(self, yy: int, mm: int, dd: int, hh: int, mt: int, ss: int):
        if not self.ValidDate(yy, mm, dd):
            return False

        if hh < 0 or hh >= 24:
            return False
        if mt < 0 or mt >= 60:
            return False
        if ss < 0 or ss >= 60:
            return False

        yp = yy + ((mm - 3) // 10)

        if (yy > 1582) or (yy == 1582 and mm > 10) or (yy == 1582 and mm == 10 and dd >= 15):
            init = 1721119.5
            jdy = yp * 365.25 - (yp // 100) + (yp // 400)
        elif (yy < 1582) or (yy == 1582 and mm < 10) or (yy == 1582 and mm == 10 and dd <= 4):
            init = 1721117.5
            jdy = yp * 365.25
        else:
            return False

        mp = (mm + 9) % 12
        jdm = mp * 30 + ((mp + 1) * 34 // 57)
        jdd = dd - 1
        jdh = (hh + (mt + (ss / 60)) / 60) / 24
        return jdy + jdm + jdd + jdh + init

    def ValidDate(self, yy, mm, dd):
        if yy < -1000 or yy > 3000:
            return False
        if mm < 1 or mm > 12:
            return False
        if yy == 1582 and mm == 10 and dd >= 5 and dd < 15:
            return False

        ndf1 = -int(yy % 4 == 0)
        ndf2 = ((yy % 400 == 0) - (yy % 100 == 0)) and (yy > 1582)
        ndf = ndf1 + ndf2
        dom = 30 + ((abs(mm - 7.5) + 0.5) % 2) - int(mm == 2) * (2 + ndf)
        if dd <= 0 or dd > dom:
            if ndf == 0 and mm == 2 and dd == 29:
                pass
            else:
                pass
            return False

        return True

    def GetSolarDays(self, yy: int, mm: int) -> int:
        if yy < -1000 or yy > 3000:
            return 0

        if mm < 1 or mm > 12:
            return 0

        ndf1 = -int(yy % 4 == 0)
        ndf2 = (int(yy % 400 == 0) - int(yy % 100 == 0)) and (yy > 1582)
        ndf = ndf1 + ndf2
        return 30 + int((abs(mm - 7.5) + 0.5) % 2) - int(mm == 2) * (2 + ndf)

    def GetLunarDays(self, yy: int, mm: int, is_leap: int) -> int:
        if yy < -1000 or yy > 3000:
            return 0
        if mm < 1 or mm > 12:
            return 0

        jdzq, jdnm, mc = self.GetZQandSMandLunarMonthCode(yy)

        leap = 0
        for j in range(1, 15):
            if mc[j] - int(mc[j]) > 0:
                leap = int(mc[j] + 0.5)
                break

        mm += 2

        nofd = [0] * 15
        for i in range(15):
            nofd[i] = int(jdnm[i + 1] + 0.5) - int(jdnm[i] + 0.5)

        dy = 0
        er = 0

        if is_leap:
            if leap < 3:
                er = 1
            else:
                if leap != mm:
                    er = 2
                else:
                    dy = nofd[mm]
        else:
            if leap == 0:
                dy = nofd[mm - 1]
            else:
                dy = nofd[mm + (mm > leap) - 1]

        return int(dy)

    def GetLeap(self, yy: int):
        jdzq, jdnm, mc = self.GetZQandSMandLunarMonthCode(yy)
        leap = 0
        for j in range(1, 15):
            if mc[j] - int(mc[j]) > 0:
                leap = int(mc[j] + 0.5)
                break
        return max(0, leap - 2)

    def Lunar2Solar(self, yy: int, mm: int, dd: int, isleap):
        if yy < -7000 or yy > 7000:
            return False
        if yy < -1000 or yy > 3000:
            return False
        if mm < 1 or mm > 12:
            return False
        if dd < 1 or dd > 30:
            return False

        jdzq, jdnm, mc = self.GetZQandSMandLunarMonthCode(yy)

        leap = 0
        for j in range(1, 15):
            if mc[j] - int(mc[j]) > 0:
                leap = int(mc[j] + 0.5)
                break

        mm = mm + 2

        nofd = []
        for i in range(15):
            nofd.append(int(jdnm[i + 1] + 0.5) - int(jdnm[i] + 0.5))

        jd = 0
        er = 0

        if isleap:
            if leap < 3:
                er = 1
            else:
                if leap != mm:
                    er = 2
                else:
                    if dd <= nofd[mm]:
                        jd = jdnm[mm] + dd - 1
                    else:
                        er = 3
        else:
            if leap == 0:
                if dd <= nofd[mm - 1]:
                    jd = jdnm[mm - 1] + dd - 1
                else:
                    er = 4
            else:
                if dd <= nofd[mm + (mm > leap) - 1]:
                    jd = jdnm[mm + (mm > leap) - 1] + dd - 1
                else:
                    er = 4

        return False if er else self.Julian2Solar(jd)[:3]

    def Solar2Lunar(self, yy: int, mm: int, dd: int):
        if not self.ValidDate(yy, mm, dd):
            return False

        prev = 0
        isleap = 0

        jdzq, jdnm, mc = self.GetZQandSMandLunarMonthCode(yy)

        jd = self.Solar2Julian(yy, mm, dd, 12, 0, 0)
        if int(jd) < int(jdnm[0] + 0.5):
            prev = 1
            jdzq, jdnm, mc = self.GetZQandSMandLunarMonthCode(yy - 1)

        for i in range(15):
            if int(jd) >= int(jdnm[i] + 0.5) and int(jd) < int(jdnm[i + 1] + 0.5):
                mi = i
                break

        if mc[mi] < 2 or prev == 1:
            yy = yy - 1

        if (mc[mi] - int(mc[mi])) * 2 + 1 != 1:
            isleap = 1

        mm = int((int(mc[mi]) + 10) % 12) + 1
        dd = int(jd) - int(jdnm[mi] + 0.5) + 1

        return [yy, mm, dd, isleap]

    def Get24JieQi(self, yy: int) -> list:
        jq = []

        dj = self.GetAdjustedJQ(yy - 1, 21, 23)
        for k, v in enumerate(dj):
            if k < 21 or k > 23:
                continue
            jq.append(self.Julian2Solar(v))

        dj = self.GetAdjustedJQ(yy, 0, 20)
        for k, v in enumerate(dj):
            jq.append(self.Julian2Solar(v))

        return jq

    def GetGanZhi(self, yy: int, mm: int, dd: int, hh: int, mt: int = 0, ss: int = 0) -> list:
        jd = self.Solar2Julian(yy, mm, dd, hh, mt, max(1, ss))
        print(jd)
        if jd is None:
            return []

        tg, dz = [], []

        jq = self.GetPureJQsinceSpring(yy)
        if jd < jq[1]:
            yy -= 1
            jq = self.GetPureJQsinceSpring(yy)

        ygz = ((yy + 4712 + 24) % 60 + 60) % 60
        tg.append(ygz % 10)
        dz.append(ygz % 12)

        for j in range(0, 15):
            if jq[j] >= jd:
                ix = j - 1
                break

        tmm = ((yy + 4712) * 12 + (ix - 1) + 60) % 60
        mgz = (tmm + 50) % 60
        tg.append(mgz % 10)
        dz.append(mgz % 12)

        jda = jd + 0.5
        thes = ((jda - int(jda)) * 86400) + 3600
        dayjd = int(jda) + thes / 86400
        dgz = (int(dayjd + 49) % 60 + 60) % 60
        tg.append(dgz % 10)
        dz.append(dgz % 12)
        if self.zwz and hh >= 23:
            tg[2] = (tg[2] + 10 - 1) % 10
            dz[2] = (dz[2] + 12 - 1) % 12

        dh = dayjd * 12
        hgz = (int(dh + 48) % 60 + 60) % 60
        tg.append(hgz % 10)
        dz.append(hgz % 12)

        return [tg, dz, jd, jq, ix]

    def GetInfo(self, gd: int, yy: int, mm: int, dd: int, hh: int, mt: int = 0, ss: int = 0) -> dict:
        if gd not in [0, 1]:
            return {}

        ret = {}
        big_tg, big_dz = [], []

        tg, dz, jd, jq, ix = self.GetGanZhi(yy, mm, dd, hh, mt, ss)

        pn = tg[0] % 2

        if (gd == 0 and pn == 0) or (gd == 1 and pn == 1):
            span = jq[ix + 1] - jd

            for i in range(1, 13):
                big_tg.append((tg[1] + i) % 10)
                big_dz.append((dz[1] + i) % 12)
        else:
            span = jd - jq[ix]

            for i in range(1, 13):
                big_tg.append((tg[1] + 20 - i) % 10)
                big_dz.append((dz[1] + 24 - i) % 12)

        days = int(span * 4 * 30)
        y = int(days / 360)
        m = int(days % 360 / 30)
        d = int(days % 360 % 30)

        ret['tg'] = tg
        ret['dz'] = dz
        ret['big_tg'] = big_tg
        ret['big_dz'] = big_dz
        ret['start_desc'] = f"{y}年{m}月{d}天起运"
        start_jdtime = jd + span * 120
        ret['start_time'] = self.Julian2Solar(start_jdtime)

        ret['bazi'] = ret['big'] = ret['years'] = ''
        ret['big_start_time'] = []

        for i in range(4):
            ret['bazi'] += self.ctg[tg[i]]
            ret['bazi'] += self.cdz[dz[i]]

        for i in range(12):
            ret['big'] += self.ctg[big_tg[i]]
            ret['big'] += self.cdz[big_dz[i]]
            ret['big_start_time'].append(self.Julian2Solar(start_jdtime + i * 10 * 360))

        j, t, d = 0, (tg[1] + 1) % 10, (dz[1] + 1) % 12
        while True:
            if (yy + j) < ret['start_time'][0]:
                j += 1
                continue
            if j >= 120:
                break

            ret['years'] += self.ctg[t]
            ret['years'] += self.cdz[d]
            if j % 10 == 0:
                ret['years'] += "\n"

            j += 1
            t = (t + 1) % 10
            d = (d + 1) % 12

        return ret

    def GetMonthsFromYear(self, ySN: int):
        ret = []
        nv = 2 + 12 * (ySN % 10)

        for i in range(12):
            pv = (i + nv) % 60
            ret.append(pv)

        return ret

    def GetHoursFromDay(self, dSN: int):
        ret = []
        nv = 2 + 12 * (dSN % 10)

        for i in range(12):
            pv = (i + nv) % 60
            ret.append(pv)

        return ret

    def ReverseBazi(self, ygz: int, mgz: int, dgz: int, hgz: int, yeai: int, mx: int):
        ret = []

        if ygz < 0 or ygz >= 60 or mgz < 0 or mgz >= 60 or dgz < 0 or dgz >= 60 or hgz < 0 or hgz >= 60:
            return ret

        if mgz not in self.GetMonthsFromYear(ygz):
            return ret

        if hgz not in self.GetHoursFromDay(dgz):
            return ret

        yeaf = yeai + mx * 60

        if yeai < -1000 or yeaf > 3000:
            return ret

        for m in (0, mx-1):
            yea = yeai + m * 60
            syc = (yea + 56) % 60
            asyc = (ygz + 60 - syc) % 60
            iy = yea + asyc
            jdpjq = self.GetPureJQsinceSpring(iy)
            mgzo = (mgz + 60 - 2) % 12
            ijd = jdpjq[mgzo]
            fjd = jdpjq[mgzo + 1]
            sdc = (math.floor(ijd) + 49) % 60
            asdc = (dgz + 60 - sdc) % 60
            idd = math.floor(ijd + asdc)
            ihh = hgz % 12
            id = idd + (ihh * 2 - 13) / 24
            fd = idd + (ihh * 2 - 11) / 24

            if fd < ijd or id > fjd:
                continue
            else:
                if id > ijd and fd < fjd:
                    ids = id
                    fds = fd

                if id < ijd and fd > ijd:
                    ids = ijd
                    fds = fd

                if id < fjd and fd > fjd:
                    ids = id
                    fds = fjd

                ret.append([self.Julian2Solar(ids), self.Julian2Solar(fds)])

        return ret
