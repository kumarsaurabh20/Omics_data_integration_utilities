"""
author: Kumar Saurabh Singh
url: https://kumarsaurabhsingh.com
creation date: 09-02-2022
"""
import sqlite3
from sqlite3 import Error

class RowsKeys:
    def __init__(self):
        self.sqlite_file = "mvc.db"
        self.conn = None            # set the placeholder for the connection
        self.create_connection()    # create the connection
        self.drop_table()           # drop the table if it exists
        self.create_table()         # creation the dummy table
        self.create_data()          # for filling up the database with dummy data

    def create_connection(self):
        """ create a database connection to the SQLite database
            specified by db_file
        :db_file: self.sqlite_file
        :creates : self.conn Connection object
        """
        try:
            self.conn = sqlite3.connect(self.sqlite_file)
            self.conn.row_factory = sqlite3.Row  #this for getting the column names!
        except Error as e:
            print("create_connection: {}".format(e))
        else:
            print("Database connection created!")

    def drop_table(self):
        """
        small function to drop the dummy table
        """
        sql = '''DROP TABLE IF EXISTS `addresses` '''
        try:
            self.conn.execute(sql)
        except Error as e:
            print("create_table: {}".format(e))
        else:
            print("Table dropped")

    def create_table(self):
        """
        small function to create a dummy table
        """
        sql = '''CREATE TABLE IF NOT EXISTS `addresses` (`id` integer PRIMARY KEY, 
                          `name` TEXT, 
                          `address` TEXT, 
                          `city` TEXT)'''
        try:
            self.conn.execute(sql)
        except Error as e:
            print("create_table: {}".format(e))
        else:
            print("Table created!")

    def create_data(self):
        addresses = [("Jansen", "Blaak 928", "Rotterdam"), ("Klaasen", "Maasberglaan 23", "Rotterdam"),
                     ("Sluijsen", "Maasstraat 25", "Barendrecht"), ("de Vos", "Meent 198", "Rotterdam"),
                     ("De Notenkraker", "Pennylane 15", "Amsterdam")]

        sql = """INSERT INTO `addresses` (`name`, `address`, `city`)
                            VALUES (?, ?, ?)"""

        try:
            cur = self.conn.cursor()
            cur.executemany(sql, addresses)
            self.conn.commit()

        except Error as e:
            print("create_table: {}".format(e))
        else:
            print("Insert of fake data!")

    def get_rows(self, fields):
        """
        Small function for getting multiple rows
        :param fields:
        :return: rows
        """
        try:
            sql = '''SELECT `name`, `address`, `city` 
                     FROM `addresses` WHERE `city` = ?'''

            cur = self.conn.cursor()
            cur.execute(sql, fields)
            return cur.fetchall()

        except Error as e:
            print("get_row: {}".format(e))

    def get_row(self, fields):
        try:
            sql = '''SELECT `name`, `address`, `city` 
                     FROM `addresses` WHERE `city` = ?'''

            cur = self.conn.cursor()
            cur.execute(sql, fields)
            return cur.fetchone()

        except Error as e:
            print("get_row: {}".format(e))

    def close_conn(self):
        try:
            self.conn.close()
        except Error as e:
            print("close_conn: {}".format(e))
        else:
            print("Connection closed!")


if __name__ == "__main__":
    s = RowsKeys()

    # get one row and print as dictionary

    print("Return one Row")

    fields = ["Barendrecht"]
    data = s.get_row(fields)
    print(dict(data))

    print("==============")

    print("Return multiple Rows")
    # get multiple rows and print as dictionary
    fields = ["Rotterdam"]
    rows = s.get_rows(fields)
    for row in rows:
        print(dict(row))

    print()
    s.close_conn()
