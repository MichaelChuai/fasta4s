/**
  * Created by michael on 1/6/16.
  */
package org.bm2.fasta4s

import java.io.{
File, RandomAccessFile, FileNotFoundException,
ObjectInputStream, ObjectOutputStream,
FileInputStream, FileOutputStream,
BufferedWriter, FileWriter
}
import java.util.NoSuchElementException
import scala.io.Source
import scala.collection.mutable.{Map=>MMap}

class Fasta {
  private var faFileName: String = ""
  private var flatFileName: String = ""
  private var gdxFileName: String = ""
  private var gdxMap: Map[String, Array[Int]] = null
  private var flatFile: RandomAccessFile = null
  def this(fileName: String) = {
    this()
    faFileName = fileName
    flatFileName = s"$faFileName.${Fasta.flatExt}"
    gdxFileName = s"$faFileName.${Fasta.gdxExt}"
    _testFile()
    _prepare()
  }

  // one based
  def apply(chr: String)(start: Int, end: Int)(implicit posStrand: Boolean = true) = {
    if (! gdxMap.contains(chr)) throw new NoSuchElementException(s"$chr not exist in this genome")
    val maxLen = gdxMap(chr)(1) - gdxMap(chr)(0)
    val demandLen = if (end >= start) end - start + 1 else 0
    val curLen = if (maxLen > demandLen) demandLen else maxLen
    val from = gdxMap(chr)(0) + start - 1
    if (posStrand) _extractSeq(from, curLen)
    else {
      _extractSeq(from, curLen).map(Fasta.ComplementCode(_)).reverse
    }
  }

  def chrom(chr: String)(implicit posStrand: Boolean = true) = {
    if (! gdxMap.contains(chr)) throw new NoSuchElementException(s"$chr not exist in this genome")
    val curLen = gdxMap(chr)(1) - gdxMap(chr)(0)
    val from = gdxMap(chr)(0)
    if (posStrand) _extractSeq(from, curLen)
    else {
      _extractSeq(from, curLen).map(Fasta.ComplementCode(_)).reverse
    }
  }

  def chromMap: Map[String, Int] = {
    gdxMap map {
      case (k, v) =>
        (k, v(1) - v(0))
    }
  }

  private def _extractSeq(from: Int, length: Int): String = {
    flatFile.seek(from)
    val c = for (i <- (1 to length).toIterator) yield flatFile.readUnsignedByte.toChar
    c.mkString
  }

  private def _prepare(): Unit = {
    val gdxFile = new ObjectInputStream(new FileInputStream(gdxFileName))
    gdxMap = gdxFile.readObject().asInstanceOf[Map[String, Array[Int]]]
    flatFile = new RandomAccessFile(flatFileName, "r")
  }

  private def _testFile(): Unit = {
    val isFlatAndGdxExist: Boolean = _testFlatAndGdx()
    val isFaExist: Boolean = _testFa()
    if (! isFlatAndGdxExist) {
      if (! isFaExist) throw new FileNotFoundException(s"$faFileName not exist")
      else _genFlatAndGdx()
    }
  }

  private def _testFa(): Boolean = {
    val faFile = new File(faFileName)
    faFile.exists
  }

  private def _testFlatAndGdx(): Boolean = {
    val flatFile = new File(flatFileName)
    val gdxFile = new File(gdxFileName)
    flatFile.exists && gdxFile.exists
  }

  private def _genFlatAndGdx(): Unit = {
    val flatFile = new BufferedWriter(new FileWriter(flatFileName))
    var gdxMap = MMap.empty[String, Array[Int]]
    var curHeader = "None"
    for (line <- Source.fromFile(faFileName).getLines if !line.trim.isEmpty) {
      line.trim match {
        case cleanLine if cleanLine.startsWith(">") =>
          val idxSep = if (curHeader == "None") 0 else gdxMap(curHeader)(1)
          gdxMap += (cleanLine.substring(1) -> Array(idxSep, idxSep))
          curHeader = cleanLine.substring(1)
        case cleanLine =>
          flatFile.write(cleanLine)
          gdxMap(curHeader)(1) += cleanLine.size
      }
    }
    flatFile.close()
    val gdxFile = new ObjectOutputStream(new FileOutputStream(gdxFileName))
    gdxFile.writeObject(gdxMap.toMap)
    gdxFile.close()
  }
}

object Fasta {
  private val flatExt = "flat"
  private val gdxExt = "gdx"
  private val ComplementCode: Map[Char, Char] = Map(
    'A' -> 'T',
    'T' -> 'A',
    'C' -> 'G',
    'G' -> 'C',
    'a' -> 't',
    't' -> 'a',
    'c' -> 'g',
    'g' -> 'c',
    'N' -> 'N',
    'n' -> 'n',
    'X' -> 'X',
    'x' -> 'x'
  )

  def apply(fileName: String) = {
    new Fasta(fileName)
  }

}
