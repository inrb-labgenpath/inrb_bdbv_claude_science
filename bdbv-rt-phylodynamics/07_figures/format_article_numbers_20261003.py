#!/usr/bin/env python3
"""format_article_numbers_20261003.py - turns a FORMAT and the VALUES of a quantity into the token that the article prints.

Standard library only. One thread. Reads no file unless it is called with 'apply'.

    python format_article_numbers_20261003.py --self-test            self-tests on made-up values (exit code 0 = all pass)
    python format_article_numbers_20261003.py --self-test LOG.csv    the same, and the table of the tests is written
    python format_article_numbers_20261003.py apply FORMATS.csv OUT.csv [COLUMN OF THE VALUES]
        FORMATS.csv holds the columns no, format and values (both JSON; the column of the values can be named);
        OUT.csv holds no, token, error

THE FORMAT (a JSON object)
    kind      a word for the kind of format (used for counting only; the program does not read it)
    pattern   the token with every value replaced by a slot <<1>>, <<2>>, ... ; all other characters of the pattern are
              printed as they stand (the dash of an interval, brackets, ' %', ' x 10^{-3}', words such as 'to')
    slots     a list; slot i fills <<i>> and takes value i of the list of values. Fields of a slot:
      part        what the slot prints (median, lower bound, upper bound, first and last day, count, ...)
      type        'number' | 'date' | 'range of dates' | 'word' | 'text'
     type 'number'
      decimals    digits after the decimal point (0 = whole number)
      factor      the value is multiplied by this before it is rounded ('1', '100' for a share printed in per cent,
                  '0.1', '1/7' for days printed in weeks, ...)
      thousands   ',' or '' (separator of the groups of three digits before the decimal point)
      sign        'minus' (a negative value carries the minus sign) | 'plus and minus' (a positive value carries '+')
                  | 'none' (the absolute value is printed; the sign is said by a word of the sentence)
      rounding    'half even' (ties go to the even digit; the digits are those of the stored text)
                  'half up'   (ties go away from zero)
                  'binary'    (the stored text is read as a double and printed as Python prints it)
                  'up' | 'down' (towards plus | minus infinity: for a bound such as 'at most 0.16' or '0.92 or more')
      limits      optional {'above': '0.99', 'below': '0.01'}: a value above (below) the limit prints the limit itself;
                  the sign or word of the bound ('>' or 'above') is not part of the token
      step        optional: the value is rounded to a multiple of this (for a bound such as '15,000 or more')
      scientific  optional true: mantissa x 10^{exponent}; 'decimals' are those of the mantissa
      transform   optional 'reciprocal': 1 / value is printed (a setting stored as a scale and printed as a rate)
      several     optional 'print once when alike': the value is a list; every member is formatted, and the token is
                  printed when all members give the same token; otherwise the format stops (the sentence prints ONE
                  value for several quantities and has to be read by a person)
      a rounded value of zero carries no sign
     type 'date'
      day         true | false ; month 'long' | 'short' | 'none' ; year true | false
     type 'range of dates' (the value is a list of two days, first and last)
      separator   the characters between the two days (' - ' with an en dash, ' to ', ' and ')
      month       'long' | 'short' ; year 'none' | 'last' (after the last day; after both days when the years differ)
      within_a_month   'merged' (the month is printed once: '4 - 19 July') | 'not attested' (the format stops)
     type 'word'
      capital     true | false ; largest: the largest number written as a word; decimals/factor/rounding as for a number
      words       optional: a mapping of special values to words
     type 'text'
      cut         'whole' | 'first N characters'

CODE POINTS (as the Word file holds them): minus sign U+2212, dash of an interval and of a range U+2013,
multiplication sign U+00D7, the space before a per cent sign and around a dash U+0020, thousands separator U+002C.
Superscripts are written ^{...} and subscripts _{...}, as in the inventory.
"""
import sys, json, csv, re
from decimal import Decimal, ROUND_HALF_EVEN, ROUND_HALF_UP, ROUND_CEILING, ROUND_FLOOR, InvalidOperation, getcontext
from fractions import Fraction
import datetime

getcontext().prec = 60
MINUS = "\u2212"; EN_DASH = "\u2013"; TIMES = "\u00d7"; SPACE = "\u0020"
MONTHS_LONG = ["January", "February", "March", "April", "May", "June", "July", "August", "September", "October", "November", "December"]
MONTHS_SHORT = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine", "ten", "eleven", "twelve"]
ROUNDINGS = {"half even": ROUND_HALF_EVEN, "half up": ROUND_HALF_UP, "up": ROUND_CEILING, "down": ROUND_FLOOR}
SLOT_RE = re.compile(r"<<(\d+)>>")


class FormatError(Exception):
    """the format cannot print the value; the message says why"""


def to_decimal(text):
    """reads a number written as text: minus sign U+2212 or '-', thousands separators ',', exponent 'e'"""
    if isinstance(text, (int, Decimal)): return Decimal(text)
    if isinstance(text, float): return Decimal(repr(text))
    s = str(text).strip().replace(MINUS, "-")
    if re.fullmatch(r"[-+]?\d{1,3}(,\d{3})+(\.\d+)?", s): s = s.replace(",", "")
    try: return Decimal(s)
    except InvalidOperation: raise FormatError("the value %r is not a number" % (text,))


def to_factor(text):
    s = str(text if text not in (None, "") else "1").strip()
    m = re.fullmatch(r"(\S+)\s*/\s*(\S+)", s)
    if m: return Decimal(m.group(1)) / Decimal(m.group(2))
    return to_decimal(s)


def to_date(text):
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})", str(text).strip())
    if not m: raise FormatError("the value %r is not a day written YYYY-MM-DD" % (text,))
    try: return datetime.date(int(m.group(1)), int(m.group(2)), int(m.group(3)))
    except ValueError: raise FormatError("the value %r is not a day of the calendar" % (text,))


def round_decimal(d, decimals, rounding="half even", stored_text=None):
    """the rounded value as a Decimal with exactly 'decimals' digits after the point"""
    q = Decimal(1).scaleb(-int(decimals))
    if rounding == "binary":
        return Decimal(format(float(d), ".%df" % int(decimals)))
    if rounding not in ROUNDINGS: raise FormatError("the rule of rounding %r is not known" % (rounding,))
    return d.quantize(q, rounding=ROUNDINGS[rounding])


def group_thousands(digits, sep):
    if not sep or len(digits) < 4: return digits
    out = []
    while len(digits) > 3: out.insert(0, digits[-3:]); digits = digits[:-3]
    out.insert(0, digits); return sep.join(out)


def print_decimal(r, decimals, thousands="", sign="minus"):
    """prints a rounded Decimal: sign, groups of three, decimal point"""
    neg = r < 0; zero = (r == 0)
    s = format(abs(r), "f")
    if "." in s: ip, fp = s.split(".")
    else: ip, fp = s, ""
    fp = (fp + "0" * int(decimals))[:int(decimals)]
    body = group_thousands(ip, thousands) + ("." + fp if int(decimals) > 0 else "")
    if zero or sign == "none": return body
    if neg: return MINUS + body
    if sign == "plus and minus": return "+" + body
    return body


def format_number(value, slot):
    if slot.get("several"):
        if not isinstance(value, (list, tuple)) or len(value) < 2: raise FormatError("the slot prints several values once; a list of two or more values is needed")
        s2 = dict(slot); s2.pop("several")
        toks = [format_number(v, s2) for v in value]
        if len(set(toks)) != 1: raise FormatError("the values do not print alike (%s); the sentence prints one value for all of them" % " | ".join(toks))
        return toks[0]
    if isinstance(value, (list, tuple)): raise FormatError("one value is needed, %d were given" % len(value))
    d = to_decimal(value)
    if slot.get("transform") == "reciprocal":
        if d == 0: raise FormatError("zero has no reciprocal")
        d = Decimal(1) / d
    elif slot.get("transform"): raise FormatError("the transform %r is not known" % (slot.get("transform"),))
    d = d * to_factor(slot.get("factor", "1"))
    decimals = int(slot.get("decimals", 0)); rounding = slot.get("rounding", "half even")
    lim = slot.get("limits") or {}
    if "above" in lim and d > to_decimal(lim["above"]): d = to_decimal(lim["above"]); rounding = "half even"
    elif "below" in lim and d < to_decimal(lim["below"]): d = to_decimal(lim["below"]); rounding = "half even"
    if slot.get("scientific"):
        if d == 0: raise FormatError("zero has no power of ten")
        e = d.adjusted(); m = round_decimal(d.scaleb(-e), decimals, rounding)
        if abs(m) >= 10: e += 1; m = round_decimal(d.scaleb(-e), decimals, rounding)
        return print_decimal(m, decimals, "", slot.get("sign", "minus")) + SPACE + TIMES + SPACE + "10^{" + (MINUS if e < 0 else "") + str(abs(e)) + "}"
    if slot.get("step"):
        st = to_decimal(slot["step"])
        if st <= 0: raise FormatError("the step must be above zero")
        d = round_decimal(d / st, 0, rounding) * st
    r = round_decimal(d, decimals, rounding)
    return print_decimal(r, decimals, slot.get("thousands", ""), slot.get("sign", "minus"))


def print_date(day, month="long", year=False, with_day=True):
    names = {"long": MONTHS_LONG, "short": MONTHS_SHORT}
    parts = []
    if with_day: parts.append(str(day.day))
    if month in names: parts.append(names[month][day.month - 1])
    elif month != "none": raise FormatError("the form of the month %r is not known" % (month,))
    if year: parts.append(str(day.year))
    if not parts: raise FormatError("the format of the date prints nothing")
    return SPACE.join(parts)


def format_date(value, slot):
    if isinstance(value, (list, tuple)): raise FormatError("one day is needed, %d were given" % len(value))
    return print_date(to_date(value), slot.get("month", "long"), bool(slot.get("year", False)), bool(slot.get("day", True)))


def format_date_range(value, slot):
    if not isinstance(value, (list, tuple)) or len(value) != 2: raise FormatError("a range of dates needs two days, first and last")
    a, b = to_date(value[0]), to_date(value[1])
    if b < a: raise FormatError("the last day lies before the first day")
    sep = slot.get("separator", SPACE + EN_DASH + SPACE); month = slot.get("month", "long"); year = slot.get("year", "none")
    if year not in ("none", "last"): raise FormatError("the form of the year %r is not known" % (year,))
    if a.year == b.year and a.month == b.month:
        how = slot.get("within_a_month", "not attested")
        if how == "merged":
            return str(a.day) + sep + print_date(b, month, year == "last")
        if how == "repeated":
            return print_date(a, month, False) + sep + print_date(b, month, year == "last")
        raise FormatError("both days fall in one month; the article holds no range of this form within one month")
    if a.year == b.year:
        return print_date(a, month, False) + sep + print_date(b, month, year == "last")
    return print_date(a, month, year == "last") + sep + print_date(b, month, year == "last")


def format_word(value, slot):
    words = slot.get("words") or {}
    if slot.get("several"):
        if not isinstance(value, (list, tuple)) or len(value) < 2: raise FormatError("the slot prints several values once; a list of two or more values is needed")
        s2 = dict(slot); s2.pop("several"); toks = [format_word(v, s2) for v in value]
        if len(set(toks)) != 1: raise FormatError("the values do not give the same word (%s)" % " | ".join(toks))
        return toks[0]
    if isinstance(value, (list, tuple)): raise FormatError("one value is needed, %d were given" % len(value))
    if str(value) in words: w = words[str(value)]
    else:
        d = to_decimal(value) * to_factor(slot.get("factor", "1"))
        r = round_decimal(d, 0, slot.get("rounding", "half even"))
        if "exact" in slot and slot["exact"] and r != d: raise FormatError("the value %s is no whole number" % value)
        n = int(r); largest = int(slot.get("largest", 10))
        if n < 0 or n > largest or n >= len(WORDS): raise FormatError("the number %d is not written as a word (the format writes words up to %d)" % (n, largest))
        w = WORDS[n]
    return (w[:1].upper() + w[1:]) if slot.get("capital") else w


def format_text(value, slot):
    if isinstance(value, (list, tuple)): raise FormatError("one text is needed, %d were given" % len(value))
    cut = slot.get("cut", "whole"); s = str(value)
    if cut == "whole": return s
    m = re.fullmatch(r"first (\d+) characters", cut)
    if m:
        n = int(m.group(1))
        if len(s) < n: raise FormatError("the text holds fewer than %d characters" % n)
        return s[:n]
    raise FormatError("the cut %r is not known" % (cut,))


FORMATTERS = {"number": format_number, "date": format_date, "range of dates": format_date_range, "word": format_word, "text": format_text}


def format_token(fmt, values):
    """fmt: the format (dict or JSON text); values: a list with one entry per slot. Returns the token."""
    if isinstance(fmt, str): fmt = json.loads(fmt)
    if isinstance(values, str): values = json.loads(values)
    slots = fmt.get("slots", []); pattern = fmt.get("pattern")
    if pattern is None: raise FormatError("the format holds no pattern")
    named = [int(x) for x in SLOT_RE.findall(pattern)]
    if sorted(named) != list(range(1, len(slots) + 1)): raise FormatError("the pattern names the slots %s, the format holds %d slots" % (named, len(slots)))
    if values is None or len(values) != len(slots): raise FormatError("the format holds %d slots, %d values were given" % (len(slots), 0 if values is None else len(values)))
    toks = {}
    for i, (sl, v) in enumerate(zip(slots, values), 1):
        if v is None or v == "" or v == []: raise FormatError("slot %d has no value" % i)
        t = sl.get("type")
        if t not in FORMATTERS: raise FormatError("the type %r of slot %d is not known" % (t, i))
        toks[i] = FORMATTERS[t](v, sl)
    return SLOT_RE.sub(lambda m: toks[int(m.group(1))], pattern)


def three_roundings(value, decimals, factor="1"):
    """the token of a number under the three common rules (half even, half up, binary); used to list the ties"""
    out = {}
    for r in ("half even", "half up", "binary"):
        out[r] = format_number(value, {"decimals": decimals, "factor": factor, "rounding": r})
    return out


# ------------------------------------------------------------------------------------------------- self-tests
def _tests():
    """every test: (name, kind of format, format, values, expected token or 'ERROR', check of code points or None).
    ALL VALUES ARE MADE UP; none is a value of the article."""
    N = lambda **k: dict(dict(type="number", decimals=2, factor="1", thousands="", sign="minus", rounding="half even"), **k)
    F1 = lambda slot, pat="<<1>>": {"kind": "test", "pattern": pat, "slots": [slot]}
    iv = {"kind": "value with interval", "pattern": "<<1>> (<<2>>" + EN_DASH + "<<3>>)", "slots": [N(), N(), N()]}
    R = lambda **k: dict(dict(type="range of dates", separator=SPACE + EN_DASH + SPACE, month="long", year="none", within_a_month="merged"), **k)
    T = []
    a = T.append
    a(("number, 2 decimals", "number with decimals", F1(N()), ["3.14159"], "3.14", None))
    a(("number, 0 decimals", "whole number", F1(N(decimals=0)), ["41.6"], "42", None))
    a(("number, 3 decimals, zeros kept", "number with decimals", F1(N(decimals=3)), ["7.1"], "7.100", None))
    # half-way values
    a(("half-way 0.125 at 2 decimals, half even", "half-way value", F1(N()), ["0.125"], "0.12", None))
    a(("half-way 0.125 at 2 decimals, half up", "half-way value", F1(N(rounding="half up")), ["0.125"], "0.13", None))
    a(("half-way 0.125 at 2 decimals, binary", "half-way value", F1(N(rounding="binary")), ["0.125"], "0.12", None))
    a(("half-way 0.135 at 2 decimals, half even", "half-way value", F1(N()), ["0.135"], "0.14", None))
    a(("half-way 2.675 at 2 decimals, half even (digits as stored)", "half-way value", F1(N()), ["2.675"], "2.68", None))
    a(("half-way 2.675 at 2 decimals, half up", "half-way value", F1(N(rounding="half up")), ["2.675"], "2.68", None))
    a(("half-way 2.675 at 2 decimals, binary (the double lies below the tie)", "half-way value", F1(N(rounding="binary")), ["2.675"], "2.67", None))
    a(("half-way 2.5 at 0 decimals, half even", "half-way value", F1(N(decimals=0)), ["2.5"], "2", None))
    a(("half-way 2.5 at 0 decimals, half up", "half-way value", F1(N(decimals=0, rounding="half up")), ["2.5"], "3", None))
    a(("half-way 3.5 at 0 decimals, half even", "half-way value", F1(N(decimals=0)), ["3.5"], "4", None))
    a(("half-way -0.125 at 2 decimals, half up (away from zero)", "half-way value", F1(N(rounding="half up")), ["-0.125"], MINUS + "0.13", None))
    # negative value that rounds to zero
    a(("negative value that rounds to zero, 2 decimals", "value that rounds to zero", F1(N()), ["-0.004"], "0.00", {"no minus sign": lambda t: MINUS not in t and "-" not in t}))
    a(("negative value that rounds to zero, with sign", "value that rounds to zero", F1(N(sign="plus and minus", decimals=3)), ["-0.0004"], "0.000", {"no sign": lambda t: t[0] == "0"}))
    a(("positive value that rounds to zero, with sign", "value that rounds to zero", F1(N(sign="plus and minus", decimals=3)), ["0.0004"], "0.000", {"no sign": lambda t: t[0] == "0"}))
    a(("negative value that rounds to zero, whole number", "value that rounds to zero", F1(N(decimals=0)), ["-0.4"], "0", None))
    # minus sign and dash
    a(("minus sign U+2212", "number with sign", F1(N(decimals=1)), ["-7.46"], MINUS + "7.5", {"first character is U+2212": lambda t: ord(t[0]) == 0x2212, "no hyphen-minus": lambda t: "-" not in t}))
    a(("minus sign U+2212 in the stored text", "number with sign", F1(N(decimals=1)), [MINUS + "7.46"], MINUS + "7.5", {"first character is U+2212": lambda t: ord(t[0]) == 0x2212}))
    a(("plus sign for a positive value", "number with sign", F1(N(decimals=3, sign="plus and minus")), ["0.0412"], "+0.041", {"first character is U+002B": lambda t: ord(t[0]) == 0x2B}))
    a(("minus sign with 'plus and minus'", "number with sign", F1(N(decimals=3, sign="plus and minus")), ["-0.0412"], MINUS + "0.041", {"first character is U+2212": lambda t: ord(t[0]) == 0x2212}))
    a(("sign said by a word of the sentence", "number with sign", F1(N(decimals=0, sign="none")), ["-33.3"], "33", None))
    a(("value with interval, dash U+2013", "value with interval", iv, ["1.2345", "0.9876", "1.5555"], "1.23 (0.99" + EN_DASH + "1.56)", {"the dash is U+2013": lambda t: t.count(EN_DASH) == 1 and "-" not in t and "\u2014" not in t, "space before the bracket is U+0020": lambda t: t[4] == SPACE}))
    a(("interval with a negative lower bound", "value with interval", {"kind": "t", "pattern": "<<1>>, 95 % HPD <<2>> to <<3>>", "slots": [N(), N(), N()]}, ["0.4", "-0.555", "1.2"], "0.40, 95 % HPD " + MINUS + "0.56 to 1.20", {"minus is U+2212": lambda t: MINUS in t and "-" not in t}))
    a(("bounds alone", "interval alone", {"kind": "t", "pattern": "<<1>>" + EN_DASH + "<<2>>", "slots": [N(), N()]}, ["0.111", "0.999"], "0.11" + EN_DASH + "1.00", None))
    # per cent
    a(("per cent from a share, space U+0020 before the sign", "per cent", F1(N(decimals=0, factor="100"), "<<1>> %"), ["0.8765"], "88 %", {"character before % is U+0020": lambda t: ord(t[-2]) == 0x20, "no no-break space": lambda t: "\u00a0" not in t and "\u202f" not in t}))
    a(("per cent, 1 decimal", "per cent", F1(N(decimals=1, factor="100"), "<<1>> %"), ["0.12345"], "12.3 %", None))
    a(("per cent stored as per cent", "per cent", F1(N(decimals=0), "<<1>> %"), ["35"], "35 %", None))
    # thousands
    a(("thousands separator, 7 digits", "thousands separator", F1(N(decimals=0, thousands=",")), ["1234567"], "1,234,567", {"separator is U+002C": lambda t: t[1] == ","}))
    a(("thousands separator, 4 digits", "thousands separator", F1(N(decimals=0, thousands=",")), ["1000"], "1,000", None))
    a(("thousands separator, 3 digits", "thousands separator", F1(N(decimals=0, thousands=",")), ["999"], "999", None))
    a(("no separator (a year)", "thousands separator", F1(N(decimals=0)), ["2031"], "2031", None))
    a(("thousands separator with decimals and sign", "thousands separator", F1(N(decimals=1, thousands=",")), ["-12345.678"], MINUS + "12,345.7", None))
    a(("stored text with separators", "thousands separator", F1(N(decimals=0, thousands=",")), ["18,765"], "18,765", None))
    a(("rounding carries into a new group", "thousands separator", F1(N(decimals=0, thousands=",")), ["999.6"], "1,000", None))
    # dates
    a(("date, long month", "date", F1(dict(type="date", month="long", year=False)), ["2031-07-04"], "4 July", None))
    a(("date, long month and year", "date", F1(dict(type="date", month="long", year=True)), ["2031-07-04"], "4 July 2031", None))
    a(("date, short month", "date", F1(dict(type="date", month="short", year=False)), ["2031-09-21"], "21 Sep", None))
    a(("date, short month and year", "date", F1(dict(type="date", month="short", year=True)), ["2031-12-01"], "1 Dec 2031", None))
    a(("date, no day of that month", "date", F1(dict(type="date", month="long", year=False)), ["2031-02-30"], "ERROR", None))
    a(("range of dates within a month", "range of dates", F1(R()), [["2031-07-04", "2031-07-19"]], "4 " + EN_DASH + " 19 July", {"dash U+2013 between two U+0020": lambda t: t[1:4] == SPACE + EN_DASH + SPACE}))
    a(("range of dates across months", "range of dates", F1(R()), [["2031-07-04", "2031-08-02"]], "4 July " + EN_DASH + " 2 August", None))
    a(("range of dates across months, year at the end", "range of dates", F1(R(year="last")), [["2031-07-04", "2031-08-02"]], "4 July " + EN_DASH + " 2 August 2031", None))
    a(("range of dates within a month, year at the end", "range of dates", F1(R(year="last")), [["2031-07-04", "2031-07-19"]], "4 " + EN_DASH + " 19 July 2031", None))
    a(("range of dates across years", "range of dates", F1(R(year="last")), [["2030-12-28", "2031-01-03"]], "28 December 2030 " + EN_DASH + " 3 January 2031", None))
    a(("range of dates, short months, across years", "range of dates", F1(R(year="last", month="short")), [["2030-12-28", "2031-01-03"]], "28 Dec 2030 " + EN_DASH + " 3 Jan 2031", None))
    a(("range of dates with 'to'", "range of dates", F1(R(separator=" to ", within_a_month="not attested")), [["2031-07-04", "2031-08-02"]], "4 July to 2 August", None))
    a(("range of dates with 'to' within a month: the format stops", "range of dates", F1(R(separator=" to ", within_a_month="not attested")), [["2031-07-04", "2031-07-19"]], "ERROR", None))
    a(("range of dates with 'and' in a sentence", "range of dates", F1(R(separator=" and ", within_a_month="not attested"), "between <<1>>"), [["2031-07-04", "2031-08-02"]], "between 4 July and 2 August", None))
    a(("range of dates, last day before first day", "range of dates", F1(R()), [["2031-08-04", "2031-07-19"]], "ERROR", None))
    a(("date with interval", "date with interval", {"kind": "t", "pattern": "<<1>> (<<2>>)", "slots": [dict(type="date", month="long", year=True), R(year="last")]}, ["2031-02-11", ["2031-01-07", "2031-03-09"]], "11 February 2031 (7 January " + EN_DASH + " 9 March 2031)", None))
    # words
    a(("number as a word", "word", F1(dict(type="word", capital=False, largest=10)), ["6"], "six", None))
    a(("number as a word, capital", "word", F1(dict(type="word", capital=True, largest=10)), ["6"], "Six", None))
    a(("number as a word, above the largest word", "word", F1(dict(type="word", capital=False, largest=10)), ["11"], "ERROR", None))
    a(("nearest whole number as a word", "word", F1(dict(type="word", capital=False, largest=10)), ["2.7"], "three", None))
    a(("days as weeks, as a word, two values alike", "word", F1(dict(type="word", capital=False, largest=10, factor="1/7", several="print once when alike")), [["20", "23"]], "three", None))
    a(("days as weeks, as a word, two values not alike", "word", F1(dict(type="word", capital=False, largest=10, factor="1/7", several="print once when alike")), [["20", "40"]], "ERROR", None))
    # bounds
    a(("bound 'above': rounded down", "bound", F1(N(rounding="down")), ["0.9987"], "0.99", None))
    a(("bound 'at most': rounded up", "bound", F1(N(rounding="up")), ["0.1512"], "0.16", None))
    a(("bound 'at most': an exact value stays", "bound", F1(N(rounding="up")), ["0.1500"], "0.15", None))
    a(("bound in per cent, rounded up", "bound", F1(N(rounding="up", decimals=0, factor="100"), "<<1>> %"), ["0.0312"], "4 %", None))
    a(("limits: value above the upper limit", "bound", F1(N(limits={"above": "0.99", "below": "0.01"})), ["0.9951"], "0.99", None))
    a(("limits: value below the lower limit", "bound", F1(N(limits={"above": "0.99", "below": "0.01"})), ["0.0003"], "0.01", None))
    a(("limits: value between the limits", "bound", F1(N(limits={"above": "0.99", "below": "0.01"})), ["0.5049"], "0.50", None))
    a(("limits: value that rounds to the limit but does not exceed it", "bound", F1(N(limits={"above": "0.99", "below": "0.01"})), ["0.9870"], "0.99", None))
    a(("bound with a step, rounded down", "bound", F1(N(rounding="down", decimals=0, step="1000", thousands=",")), ["23456.7"], "23,000", None))
    # powers of ten
    a(("unit 10^-3 in the pattern", "number, unit 10^-3", F1(N(factor="1000"), "<<1>> " + TIMES + " 10^{" + MINUS + "3}"), ["0.001234"], "1.23 " + TIMES + " 10^{" + MINUS + "3}", {"multiplication sign U+00D7": lambda t: TIMES in t and "x" not in t, "minus U+2212 in the exponent": lambda t: MINUS in t and "-" not in t}))
    a(("value with interval in units of 10^-3", "value with interval, unit 10^-3", {"kind": "t", "pattern": "<<1>> (<<2>>" + EN_DASH + "<<3>>) " + TIMES + " 10^{" + MINUS + "3}", "slots": [N(), N(), N()]}, ["1.2345", "1.1", "1.4"], "1.23 (1.10" + EN_DASH + "1.40) " + TIMES + " 10^{" + MINUS + "3}", None))
    a(("scientific notation, 1 decimal", "power of ten", F1(N(decimals=1, scientific=True)), ["250000000"], "2.5 " + TIMES + " 10^{8}", None))
    a(("scientific notation, whole mantissa", "power of ten", F1(N(decimals=0, scientific=True)), ["4,000,000"], "4 " + TIMES + " 10^{6}", None))
    a(("scientific notation, negative exponent", "power of ten", F1(N(decimals=1, scientific=True)), ["0.00025"], "2.5 " + TIMES + " 10^{" + MINUS + "4}", None))
    a(("scientific notation, mantissa rounds to ten", "power of ten", F1(N(decimals=1, scientific=True)), ["9.96e5"], "1.0 " + TIMES + " 10^{6}", None))
    a(("factor 0.1 (stored in units of 10^-4, printed in units of 10^-3)", "number, unit 10^-3", F1(N(factor="0.1")), ["6.4"], "0.64", None))
    a(("factor 0.000001 (printed in millions)", "whole number", F1(N(decimals=0, factor="0.000001"), "<<1>> million"), ["25000000"], "25 million", None))
    a(("reciprocal of the stored value", "number with decimals", F1(N(decimals=3, transform="reciprocal")), ["250"], "0.004", None))
    a(("reciprocal of zero", "errors of use", F1(N(decimals=3, transform="reciprocal")), ["0"], "ERROR", None))
    # several values printed once
    a(("two values that print alike", "one token for several values", F1(N(several="print once when alike")), [["0.7012", "0.6991"]], "0.70", None))
    a(("two values that do not print alike: the format stops", "one token for several values", F1(N(several="print once when alike")), [["0.7012", "0.7191"]], "ERROR", None))
    # range of values
    a(("range of values with 'to'", "range of values", {"kind": "t", "pattern": "<<1>> to <<2>>", "slots": [N(decimals=1), N(decimals=1)]}, ["1.44", "2.26"], "1.4 to 2.3", None))
    a(("range of values with 'between ... and'", "range of values", {"kind": "t", "pattern": "between <<1>> and <<2>>", "slots": [N(), N()]}, ["1.444", "2.266"], "between 1.44 and 2.27", None))
    # text
    a(("text, first 8 characters", "text", F1(dict(type="text", cut="first 8 characters")), ["abcdef0123456789"], "abcdef01", None))
    a(("text, whole", "text", F1(dict(type="text", cut="whole")), ["9.8.7"], "9.8.7", None))
    a(("text shorter than the cut", "text", F1(dict(type="text", cut="first 8 characters")), ["abc"], "ERROR", None))
    a(("fixed text (a pattern without slot)", "fixed text", {"kind": "fixed text", "pattern": "ABC_0123", "slots": []}, [], "ABC_0123", None))
    a(("fixed text given a value: the format stops", "fixed text", {"kind": "fixed text", "pattern": "ABC_0123", "slots": []}, ["1"], "ERROR", None))
    a(("reference: venue, year, volume, pages", "reference", {"kind": "reference", "pattern": "<<1>> <<2>>;<<3>>:<<4>>", "slots": [dict(type="text", cut="whole")] * 4}, ["Journal of Made-up Results", "2031", "7", "11-19"], "Journal of Made-up Results 2031;7:11-19", {"hyphen-minus U+002D in the pages, as given": lambda t: t.count("-") == 2}))
    a(("reference: venue and year", "reference", {"kind": "reference", "pattern": "<<1>> <<2>>", "slots": [dict(type="text", cut="whole")] * 2}, ["Made-up Letters", "2031"], "Made-up Letters 2031", None))
    # errors of use
    a(("a slot without value", "errors of use", F1(N()), [""], "ERROR", None))
    a(("fewer values than slots", "errors of use", iv, ["1.0", "2.0"], "ERROR", None))
    a(("a value that is no number", "errors of use", F1(N()), ["not evaluated"], "ERROR", None))
    a(("a pattern that names a slot the format does not hold", "errors of use", {"kind": "t", "pattern": "<<1>> <<2>>", "slots": [N()]}, ["1.0"], "ERROR", None))
    a(("a rule of rounding that is not known", "errors of use", F1(N(rounding="nearest")), ["1.0"], "ERROR", None))
    a(("format and values given as JSON text", "errors of use", json.dumps(F1(N())), json.dumps(["2.345"]), "2.34", None))
    return T


def kinds_tested():
    """the kinds of format that the self-tests cover"""
    return sorted(set(t[1] for t in _tests()))


def self_test(log_path=None):
    rows = []; n_fail = 0
    for name, kind, fmt, values, expected, checks in _tests():
        try: got = format_token(fmt, values); err = ""
        except FormatError as e: got = "ERROR"; err = str(e)
        ok = (got == expected); cp_ok = ""
        if ok and checks and got != "ERROR":
            bad = [k for k, f in checks.items() if not f(got)]
            cp_ok = "yes" if not bad else "no: " + "; ".join(bad)
            ok = ok and not bad
        n_fail += (not ok)
        rows.append(dict(test=name, kind=kind, format=json.dumps(fmt, ensure_ascii=False) if not isinstance(fmt, str) else fmt, values=json.dumps(values, ensure_ascii=False) if not isinstance(values, str) else values,
                         expected=expected, got=got, code_points_of_the_token=" ".join("U+%04X" % ord(c) for c in got) if got != "ERROR" else "", checks_of_code_points=cp_ok, message=err, passed="yes" if ok else "no"))
    if log_path:
        with open(log_path, "w", newline="", encoding="utf8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    return rows, n_fail


def apply_file(src, out, values_column="values"):
    n = 0; ne = 0
    with open(src, newline="", encoding="utf8") as f, open(out, "w", newline="", encoding="utf8") as g:
        w = csv.DictWriter(g, fieldnames=["no", "token", "error"], lineterminator="\n"); w.writeheader()
        for r in csv.DictReader(f):
            try:
                if r[values_column] in ("", "none"): raise FormatError("the row holds no values")
                w.writerow(dict(no=r["no"], token=format_token(r["format"], r[values_column]), error=""))
            except (FormatError, ValueError, KeyError) as e: w.writerow(dict(no=r["no"], token="", error=str(e))); ne += 1
            n += 1
    return n, ne


if __name__ == "__main__":
    a = sys.argv[1:]
    if a and a[0] == "--self-test":
        rows, n_fail = self_test(a[1] if len(a) > 1 else None)
        kinds = {}
        for r in rows: kinds.setdefault(r["kind"], [0, 0]); kinds[r["kind"]][0] += 1; kinds[r["kind"]][1] += (r["passed"] == "yes")
        for k, (n, p) in kinds.items(): print("%-22s %2d tests, %2d pass" % (k, n, p))
        print("self-tests: %d, pass: %d, fail: %d" % (len(rows), len(rows) - n_fail, n_fail))
        for r in rows:
            if r["passed"] != "yes": print("FAIL", r["test"], "| expected", repr(r["expected"]), "| got", repr(r["got"]), r["message"], r["checks_of_code_points"])
        sys.exit(1 if n_fail else 0)
    elif len(a) in (3, 4) and a[0] == "apply":
        n, ne = apply_file(a[1], a[2], *(a[3:4])); print("rows: %d, tokens: %d, errors: %d" % (n, n - ne, ne)); sys.exit(0)
    else:
        print(__doc__); sys.exit(2)
